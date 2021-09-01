#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    pdb_tools.py
#
#    Utility functions retrieving data from PDB/mmCIF files or
#    from the PDB/PDBe database.
#
#    Copyright (C) 2021 Valentina Sora 
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



# Standard library
import json
import logging as log
import requests
from tempfile import NamedTemporaryFile
# Third-party packages
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
# SLiMfast
from SLiMfast.defaults import STRUCTURE_METADATA



# Get the module logger
logger = log.getLogger(__name__)



def get_protein_segments(pdb_id, uniprot_id):
    """Get which portions (segments) of a protein, identified by its
    UniProt ID, are covered in a given PDB file, given its PDB ID.
    """

    # URL where the mapping between PDB entities in a PDB entry
    # and their corresponding UniProt data is stored on PDBe
    mapping_url = \
        "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{:s}"

    # Get the data from PDBe
    response = requests.get(mapping_url.format(pdb_id))
    
    # Convert to text and read the resulting dictionary through JSON
    response_text = json.loads(response.text)

    # Get the data about the mappings
    data = \
        response_text[pdb_id.lower()]["UniProt"][uniprot_id]["mappings"]
    
    # Create an empty dictionary to store the mapping between
    # the PDB ID and the chains corresponding to the protein
    # identified by the UniProt ID
    mappings = {}

    # For each mapping found
    for m in data:

        # Get the chain ID
        pdb_chain = m["chain_id"]

        # Get the start and end of the chain as stored in the PDB
        pdb_start = m["start"]["author_residue_number"]
        pdb_end = m["end"]["author_residue_number"]

        # Get the start and end of the chain as they would be
        # in UniProt numbering
        uniprot_start = m["unp_start"]
        uniprot_end = m["unp_end"]

        if pdb_start == "null" or pdb_end == "null":
            warnstr = \
                f"Mapping for UniProt ID {uniprot_id} residues " \
                f"{unp_start}-{unp_end} not available in PDB " \
                f"{pdb_id} chain {pdb_chain}. Therefore, this " \
                f"segment will not be considered."
            logger.warning(warnstr)
            continue

        # Create a tuple of tuples to store those ranges
        ranges = [[pdb_start, pdb_end], [uniprot_start, uniprot_end]]

        # If the chain ID was already encountered (i.e. the PDB
        # chain corresponds to multiple discontinuous protein
        # segments)
        if pdb_chain in mappings.keys():
            # Add the ranges to the mapping
            mappings[pdb_chain].append(ranges)

        # If it is the first time that the chain ID is encountered
        else:
            # Create a new list and add the mapping as first element
            mappings[pdb_chain] = [ranges]

    # Return the mappings
    return mappings


def get_pdb_data(pdb_id, model = None):
    """Get the content of a PDB file given the PDB ID.
    """
    
    # Generic URL where PDB files for single entries
    # are found
    pdb_url = "https://files.rcsb.org/view/{:s}.pdb"
    
    # Get the data from the Protein Data Bank
    response = requests.get(pdb_url.format(pdb_id))
    # If something went wrong getting the data
    if not response.ok:
        # Raise the corresponding error
        response.raise_for_status()

    # Create an empty list with the lines to be kept
    data = []

    # Set a boolean switch to signal whether a specific
    # line of the file should be saved or not
    save_line = True

    # Previous line
    prev_line = ""

    # Initialize the model to 1 (since in case there is
    # no MODEL field, all lines of the only model present
    # will be saved anyway)
    model = 1 if model is None else model

    for line in response.text.split("\n"):

        # If you encounter the line specifying the model
        if line.startswith("MODEL"):
            # Get the model number
            _, num = [i for i in line.split(" ") if i != ""]
            # If it is the model of interest, save the line
            # corresponding to info about it
            if int(num) == model:
                save_line = True
            # Otherwise, signal that the line should not
            # be saved
            else:
                save_line = False

        # If you reached the end of the ATOM section
        elif line.startswith("ENDMDL"):
            # Start saving the lines again
            save_line = True

        # If the save-line switch is on, save the current line
        if save_line:
            prev_line = line
            if prev_line.startswith("ENDMDL"):
                continue
            if line != "":
                data.append(line)

    # Return the data retrieved
    return data


def get_missing_residues(pdb_data):
    """Get the missing residues in a PDB file parsing
    the REMARK 465 records.
    """

    # Create an empty dictionary to store the missing residues
    missing_residues = {}
    # Create an empty list to store the strings represening the
    # missing residues (for logging purposes)
    missing_residues_list = []

    # Boolean switch to signal when to start parsing the
    # missing residues' section
    parse = False

    for line in pdb_data:

        # If the line is a REMARK 465 record
        if line.startswith("REMARK 465"):

            # Split the line and get only non-empty items
            line = [item for item in line.split(" ") if item != ""]

            # If the line has less than 3 items, skip it
            if len(line) < 3:
                continue

            # If the third item of the line is RES
            if line[2] == "RES" or line[2] == "M":
                
                # Turn on the switch for parsing the missing
                # residues' section 
                parse = True
                # Skip the line since there are no data on it
                continue

            # If the section has started
            if parse:
                
                # Get info from the line
                remark, nr, resname, chain, resnum = line

                # Add the missing residue to the chain it
                # should have belonged to
                if chain in missing_residues.keys():
                    missing_residues[chain].append(int(resnum))
                else:
                    missing_residues[chain] = [int(resnum)]

                # Add the missing residue to the list of
                # missing residues
                missing_residues_list.append(\
                    f"{chain}-{resname}{resnum}")

    # Inform the user about the missing residues
    infostr = \
        f"The following residues are missing: " \
        f"{', '.join(missing_residues_list)}."
    logger.info(infostr)

    # Return the dictionary of missing residues
    return missing_residues


def get_structure_metadata(pdb_id):
    """Get metadata associated to a structure in the PDB 
    from its mmCIF file.
    """

    # Metadata are collected from mmCIF files and not
    # directly from the PDB files since the mmCIF format
    # is more consistent

    # Generic URL where mmCIF files for single entries
    # are found
    mmcif_url = "https://files.rcsb.org/view/{:s}.cif"
    
    # Create a dictionary to store the metadata
    metadata = {}

    # Use a temporary file to store the mmCIF file
    with NamedTemporaryFile("w") as tmp:
        
        # Get the mmCIF file
        response = requests.get(mmcif_url.format(pdb_id))
        # Write the mmCIF file to the temporary file
        tmp.write(response.text)
         
        # Create a dictionary from the header data
        mmcif_dict = MMCIF2Dict(tmp.name)

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
                    else:
                        data = None
                
                # If the data type is tuple   
                else:
                    # Join the pieces of data in a string,
                    # separated by a semicolon
                    data = \
                        ":".join(datatype(\
                            [item for item in data if item != "?"]))        
            else:
                # The field is not present
                data = None

            # Append data to the list
            metadata[col] = data

        # Return the dictionary of metadata
        return metadata