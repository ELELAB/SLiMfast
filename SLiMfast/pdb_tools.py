#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    pdb_tools.py
#
#    Utilities to retrieve data from PDB/mmCIF files or from the
#    PDB/PDBe databases.
#
#    Copyright (C) 2024 Valentina Sora 
#                       <sora.valentina1@gmail.com>
#                       Matteo Tiberti 
#                       <matteo.tiberti@gmail.com> 
#                       Elena Papaleo
#                       <elenap@cancer.dk>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program. 
#    If not, see <http://www.gnu.org/licenses/>.


#######################################################################


# Set the module's description
__doc__ = \
    "Utilities to retrieve data from PDB/mmCIF files or from the " \
    "PDB/PDBe databases."


#######################################################################


# Standard library
import json
import logging as log
import requests
from tempfile import NamedTemporaryFile
# Third-party packages
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
# SLiMfast
from . import defaults


#######################################################################


# Get the module's logger
logger = log.getLogger(__name__)


#######################################################################


# The metadata (taken from MMCIF files) retrieved for PDB IDs
# of interest
STRUCTURE_METADATA = [
    # Method used for structure determination
    ("_exptl.method",
     str, 
     defaults.DF_COLUMNS["method"]),
    # Resolution
    ("_refine.ls_d_res_high",
     float,
     defaults.DF_COLUMNS["resolution"]),
    # Deposition date
    ("_pdbx_database_status.recvd_initial_deposition_date",
     str,
     defaults.DF_COLUMNS["deposition_date"]),
    ]


#######################################################################


def get_protein_segments(struct_id,
                         uniprot_id):
    """Given a protein's UniProt ID and a structure's PDB ID, identify
    which portions (segments) of the protein are covered by the
    structure.
    
    Parameters
    ----------
    struct_id : ``str``
        The structure's PDB ID.

    uniprot_id : ``str``
        The protein's UniProt ID.
    """

    # The generic URL where the mapping between entities in a PDB
    # entry and their corresponding UniProt data are stored on PDBe
    mapping_url = \
        "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{:s}"

    #-----------------------------------------------------------------#

    # Get the data from PDBe
    response = requests.get(mapping_url.format(struct_id))

    #-----------------------------------------------------------------#

    # If something went wrong
    if not response.ok:
        if response.status_code == 404:
            errstr = \
                f"The PDB structure with ID '{struct_id}' was not found on PDBe. " \
                f"Please check if the ID is correct or if the structure is available."
        else:        # Inform the user that something went wrong
            errstr = \
                "It was not possible to retrieve the mapping between " \
                f"the PDB structure with ID '{struct_id}' and the " \
                "protein segments covered by it from PDBe. " \
                f"Error code: {response.status_code}."

            # Raise an error
            raise ValueError(errstr)

    #-----------------------------------------------------------------#

    # Convert to text and read the resulting dictionary through JSON
    response_text = json.loads(response.text)


    #-----------------------------------------------------------------#

    # Check if the struct_id and uniprot_id exist in the response
    if struct_id.lower() not in response_text or \
       "UniProt" not in response_text[struct_id.lower()] or \
       uniprot_id not in response_text[struct_id.lower()]["UniProt"]:
        print(f"Mapping data not found for UniProt ID '{uniprot_id}' and PDB ID '{struct_id}'.")
        return {}

    #------------------------------------------------------------------#

    # Get the data about the mappings
    data = \
        response_text[struct_id.lower()][\
            "UniProt"][uniprot_id]["mappings"]


    #-----------------------------------------------------------------#
        
    # Create an empty dictionary to store the mapping between the PDB
    # ID and the chains corresponding to the protein of interest
    mappings = {}

    #-----------------------------------------------------------------#

    # For each mapping found
    for m in data:

        # Get the chain ID
        struct_chain = m["chain_id"]

        # Get the start of the chain in the structure's own numbering
        struct_start = m["start"]["residue_number"]

        # Get the start of the chain in the structure's own numbering
        struct_end = m["end"]["residue_number"]

        # Get the start of the chain in UniProt numbering
        uniprot_start = m["unp_start"]

        # Get the end of the chain in UniProt numbering
        uniprot_end = m["unp_end"]

        #-------------------------------------------------------------#

        # If no mapping was available
        if struct_start == "null" or struct_end == "null":

            # Warn the user
            warnstr = \
                "There is no mapping available for residues " \
                f"{unp_start}-{unp_end} of the protein with UniProt " \
                f"ID '{uniprot_id}' in chain '{struct_chain}' " \
                f"of the PDB structure with ID '{struct_id}'. " \
                "Therefore, this protein segment will be ignored."
            logger.warning(warnstr)

            # Go to the next mapping
            continue

        #-------------------------------------------------------------#

        # Create a tuple of tuples to store the mapping's ranges
        ranges = \
            [[struct_start, struct_end], \
             [uniprot_start, uniprot_end]]

        #-------------------------------------------------------------#

        # If the chain ID was already encountered (i.e., the PDB
        # chain corresponds to multiple discontinuous protein
        # segments)
        if struct_chain in mappings.keys():
                
            # Add the ranges to the mapping
            mappings[struct_chain].append(ranges)

        # If it is the first time that the chain ID is encountered
        else:
                
            # Create a new list and add the mapping as the list's
            # first element
            mappings[struct_chain] = [ranges]

        #-----------------------------------------------------------------#

        # Return the mappings found
        return mappings


def get_pdb_data(struct_id,
                 model = None):
    """Given a PDB structure's ID, retrieve the data associated
    with the corresponding PDB file and return them as a list.

    Parameters
    ----------
    struct_id : ``str``
        The structure's PDB ID.

    model : ``int``, optional
        The model of interest.

        If not provided, data for all models will be retrieved
        and returned.

    Returns
    -------
    struct_data : ``list``
        A list containing the data parsed from the structure's
        PDB file.
    """
    
    # The generic URL where PDB files for single entries are found
    struct_url = "https://files.rcsb.org/view/{:s}.pdb"

    #-----------------------------------------------------------------#
    
    # Get the data from the Protein Data Bank
    response = requests.get(struct_url.format(struct_id))

    #-----------------------------------------------------------------#
    
    # If something went wrong
    if not response.ok:
        
        # Warn the user
        errstr = \
            "It was not possible to retrieve the data for the " \
            f"structure with ID '{struct_id}' from the Protein Data " \
            f"Bank. Error code: {response.status_code}."
        logger.exception(errstr)

        # Raise an error
        raise Exception(errstr)

    #-----------------------------------------------------------------#

    # Create an empty list with the lines to be kept
    data = []

    #-----------------------------------------------------------------#

    # Set a boolean switch to signal whether a specific line of the
    # file should be saved or not
    save_line = True

    #-----------------------------------------------------------------#

    # Set a variable storing the content of the previous line and
    # initialize it to an empty string
    prev_line = ""

    #-----------------------------------------------------------------#

    # For each line in the data
    for line in response.text.split("\n"):

        # If no model was specified, all of them will be saved
        if model is None:

            # If the line is not empty
            if line != "":

                # Save the line
                data.append(line)

            # Continue
            continue

        #-------------------------------------------------------------#

        # Otherwise
        else:

            # If we encounter a line specifying the model
            if line.startswith("MODEL"):
                
                # Get the model's number
                _, num = [i for i in line.split(" ") if i != ""]
                
                # If it is the model of interest
                if int(num) == model:

                    # Turn the 'save line' switch on
                    save_line = True
                
                # Otherwise
                else:

                    # Turn the 'save line' switch off
                    save_line = False

            #---------------------------------------------------------#

            # If we reached the end of the 'ATOM' section
            elif line.startswith("ENDMDL"):
                
                # Turn the 'save line' switch on
                save_line = True

            #---------------------------------------------------------#

            # If the save-line switch is on
            if save_line:

                # Set the 'previous line' to the current line
                prev_line = line

                # If the previous line is the end of a model
                if prev_line.startswith("ENDMDL"):

                    # Continue
                    continue

                # If the line is not empty
                if line != "":

                    # Add it to the data
                    data.append(line)

    #-----------------------------------------------------------------#

    # Return the data
    return data


def get_missing_residues(pdb_data):
    """Get the missing residues in a PDB file by parsing the file's
    'REMARK 465' records.

    Parameters
    ----------
    pdb_data : ``list``
        A list containing the lines parsed from a PDB file.

    Returns
    -------
    missing_residues : ``dict``
        A dictionary mapping each chain in the structure to its
        missing residues.
    """

    # Create an empty dictionary to store the missing residues
    missing_residues = {}

    #-----------------------------------------------------------------#

    # Set a boolean switch to signal when to start parsing
    parse = False

    #-----------------------------------------------------------------#

    # For each line in the data
    for line in pdb_data:

        # If the line is a REMARK 465 record
        if line.startswith("REMARK 465"):

            # Split the line and get only non-empty items
            line = [item for item in line.split(" ") if item != ""]

            # If the line has less than 3 items
            if len(line) < 3:

                # Skip it
                continue

            #---------------------------------------------------------#

            # If the third item of the line is 'RES' or 'M'
            if line[2] == "RES" or line[2] == "M":
                
                # Turn on the switch for parsing
                parse = True
                
                # Skip the line since there are no data on it
                continue

            #---------------------------------------------------------#

            # If we are parsing
            if parse:
                
                # Get info from the line
                remark, nr, resname, chain, resnum = line

                # If we already have missing residues for the chain
                # the residue belongs to
                if chain in missing_residues.keys():

                    # Add the residue to the record for the chain in
                    # the final dictionary
                    missing_residues[chain].append(int(resnum))

                # Otherwise
                else:

                    # Create a new record for the chain in the final
                    # dictionary and add the residue
                    missing_residues[chain] = [int(resnum)]

    #-----------------------------------------------------------------#

    # Return the dictionary of missing residues
    return missing_residues


def get_structure_metadata(struct_id):
    """Get the metadata associated with a structure in the Protein Data
    Bank from the structure's mmCIF file.

    Parameters
    ----------
    struct_id : ``str``
        The structure's PDB ID.
    """

    # The metadata are collected from mmCIF files and not from the
    # PDB files since the mmCIF format is more consistent

    # The generic URL where mmCIF files for single entries are found
    mmcif_url = "https://files.rcsb.org/view/{:s}.cif"
    
    #-----------------------------------------------------------------#

    # Create a dictionary to store the metadata
    metadata = {}

    #-----------------------------------------------------------------#

    # Create a temporary file to store the mmCIF file
    with NamedTemporaryFile("w") as tmp:
        
        # Get the mmCIF file
        response = requests.get(mmcif_url.format(struct_id))

        #-------------------------------------------------------------#
        
        # Write the mmCIF file to the temporary file
        tmp.write(response.text)

        #-------------------------------------------------------------#
         
        # Create a dictionary from the header data
        mmcif_dict = MMCIF2Dict(tmp.name)

        #-------------------------------------------------------------#

        # For each field of interest
        for field, datatype, col in STRUCTURE_METADATA:
            
            # If the field is present in the mmCIF file
            if field in mmcif_dict.keys():

                # Retrieve data for that field
                data = mmcif_dict[field]
                
                # If the data type is not tuple
                if datatype is not tuple:
                    
                    # If the piece of data is available
                    if data[0] != "?":
                        
                        # Convert it to the appropriate data type
                        data = str(datatype(data[0]))
                    
                    # Otherwise
                    else:

                        # The data for that field will be None
                        data = None
                
                # If the data type is tuple   
                else:
                    
                    # Join the pieces of data in a string,
                    # separated by a semicolon
                    data = \
                        ":".join(datatype(\
                            [item for item in data if item != "?"]))

            #---------------------------------------------------------#   
            
            # Otherwise
            else:
                
                # The field is not present
                data = None

            #---------------------------------------------------------#  

            # Add the data to the dictionary
            metadata[col] = data

        #-------------------------------------------------------------#

        # Return the dictionary of metadata
        return metadata
