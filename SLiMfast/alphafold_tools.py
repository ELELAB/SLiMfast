#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    alphafold_tools.py
#
#    Utility functions to retrieve data from PDB/mmCIF files from
#    the AlphaFold Protein Structure Database.
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
    "Utility functions to retrieve data from PDB/mmCIF files from " \
    "the AlphaFold Protein Structure Database."


#######################################################################


# Standard library
import itertools
import json
import logging as log
import requests
from tempfile import NamedTemporaryFile
# Third-party packages
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
# SLiMfast
from . import defaults


#######################################################################


# Get the module logger
logger = log.getLogger(__name__)


#######################################################################


# The metadata (taken from MMCIF files) retrieved for the AlphaFold
# IDs of interest
STRUCTURE_METADATA = [
    
    # The beginning of the sequence covered by the structure in
    # UniProt numbering
    ("_ma_target_ref_db_details.seq_db_align_begin",
     int,
     "uniprot_start"),
    
    # The end of the sequence covered by the structure in UniProt
    # numbering
    ("_ma_target_ref_db_details.seq_db_align_end",
     int,
     "uniprot_end"),
    
    # The beginning of the sequence covered by the structure in
    # in the structure's own numbering
    ("_struct_ref.pdbx_align_begin",
     int,
     "struct_start"),

    # The end of the sequence covered by the structure in the
    # structure's own numbering
    ("_struct_ref.pdbx_align_begin",
     int,
     "struct_end"),

    # The alphaFold ID
    ("_database_2.database_code",
     str, 
     defaults.DF_COLUMNS["alphafold_id"]),

    # The chain ID
    ("_entity_poly.pdbx_strand_id",
     str,
     defaults.DF_COLUMNS["struct_chain"]),

    # The UniProt isoform
    ("_ma_target_ref_db_details.seq_db_align_isoform",
     str,
     defaults.DF_COLUMNS["uniprot_isoform"]),

    # The deposition date
    ("_pdbx_database_status.recvd_initial_deposition_date",
     str,
     defaults.DF_COLUMNS["deposition_date"]),
    ]


#######################################################################


def get_alphafold_data(uniprot_id = None,
                       struct_id = None):
    """Get the PDB structure(s) associated with a specific UniProt ID
    in the AlphaFold Protein Structure Database, or a specific
    structure from the database.

    Parameters
    ----------
    uniprot_id : ``str``, optional
        The UniProt ID of the protein whose structure(s) should be
        retrieved from the AlphaFold Protein Structure Database.

        It needs to be provided if ``struct_id`` is not provided.

    struct_id : ``str``, optional
        The ID of a specific structure stored in the AlphaFold
        Protein Structure Database.

        It needs to be provided if ``uniprot_id`` is not provided.

    Returns
    -------
    alphafold_data : ``dict``
        A dictionary mapping each structure's ID to a tuple
        containing:

        - The structure's version.
        - The structure's fragment.
        - The content of the PDB file associated with the
          structure.
    """

    # If neither the UniProt ID nor the structure ID were passed
    if uniprot_id is None and struct_id is None:

        # Raise an error
        errstr = \
            "You must provide either a UniProt ID or a specific " \
            "structure's ID (i.e., AF-P10415-F1-v3)."
        raise ValueError(errstr)

    #-----------------------------------------------------------------#

    # Set the generic string representing the ID of the structures
    # deposited in AlphaFold Protein Structure Database
    struct_id_fmt = "AF-{:s}-F{:d}-model_v{:d}"

    # Set the generic URL where AlphaFoldDB PDB files are found
    struct_url = "https://www.alphafold.ebi.ac.uk/files/{:s}.pdb"

    #-----------------------------------------------------------------#

    # If the structure ID was specified
    if struct_id is not None:

        # Try to get the data from the AlphaFold Protein Structure
        # Database
        response = requests.get(struct_url.format(struct_id))

        # If something went wrong
        if not response.ok:

            # Re-raise the exception
            response.raise_for_status()

        # Return the data
        return [l for l in response.text.split("\n")]

    #-----------------------------------------------------------------#
    
    # Create an empty dictionary to store the structures
    alphafold_data = {}

    #-----------------------------------------------------------------#

    # Try to get all versions of the structure
    for version in itertools.count(start = 6):

        # Try to get all fragments (if the structure is split into
        # multiple fragments)
        for fragment in itertools.count(start = 1):

            # Set the URL where to retrieve the structure
            curr_struct_url = \
                struct_url.format(struct_id_fmt.format(uniprot_id,
                                                       fragment,
                                                       version))
            
            # Try to get the data from the AlphaFold Protein
            # Structure Database
            response = requests.get(curr_struct_url)

            #---------------------------------------------------------#

            # If the structure was not found (it does
            # not exist)
            if response.status_code == 404:

                # Break out of the loop since we do not have more
                # fragments
                break

            #---------------------------------------------------------#

            # If the structure was found
            else:

                # If everything went smoothly
                if response.ok:

                    # Get the ID of the structure
                    struct_id = \
                        struct_id_fmt.format(uniprot_id,
                                             fragment,
                                             version)

                    # Add the data to the dictionary
                    alphafold_data[struct_id] = \
                        (version,
                         fragment, \
                         [l for l in response.text.split("\n")])

                # If something else went wrong
                else:

                    # Raise the corresponding error
                    response.raise_for_status()

        #-------------------------------------------------------------#

        # If the structure was not found and it is because we reached
        # he last version available
        if response.status_code == 404 and fragment == 1:

            # Break out of the loop
            break

    #-----------------------------------------------------------------#

    # Return the data retrieved
    return alphafold_data


def get_structure_metadata(struct_id):
    """Get metadata associated with a structure in the AlphaFold
    Protein Structure Database from the structure's mmCIF file.
    """

    # The metadata are collected from mmCIF files and not
    # directly from the PDB files since the mmCIF format
    # is more recent, more consistent, and contains
    # more information

    # Set the generic URL where mmCIF files for single entries
    # are found
    mmcif_url = "https://www.alphafold.ebi.ac.uk/files/{:s}.cif"

    #-----------------------------------------------------------------#
    
    # Create a dictionary to store the metadata
    metadata = {}

    #-----------------------------------------------------------------#

    # Open a temporary file to store the mmCIF file
    with NamedTemporaryFile("w") as tmp:
        
        # Try to get the mmCIF file
        response = requests.get(mmcif_url.format(struct_id))

        # If something went wrong
        if not response.ok:

            # Raise an error
            errstr = \
                "It was not possible to retrieve the mmCIF " \
                f"file for the structure with ID '{struct_id}' " \
                "from the AlphaFold Protein Structure Database " \
                "(to get the structure's metadata). Error code: " \
                f"{response.status_code}."
            raise Exception(errstr)

        #-------------------------------------------------------------#

        # Write the mmCIF file to the temporary file
        tmp.write(response.text)

        #-------------------------------------------------------------#
         
        # Create a dictionary from the mmCIF file's header
        mmcif_dict = MMCIF2Dict(tmp.name)

        #-------------------------------------------------------------#

        # Initialize a list to store the segment
        # covered by the AlphaFold structure in the
        # UniProt numbering
        unp_segment = [None, None]

        # Initialize a list to store the segment
        # covered by the AlphaFold structure in the PDB
        # numbering
        struct_segment = [None, None]

        #-------------------------------------------------------------#

        # For each field of interest
        for field, datatype, col in STRUCTURE_METADATA:

            # If the field is present in the mmCIF file
            if field in mmcif_dict.keys():

                # Retrieve the data for that field
                data = mmcif_dict[field]

                # If the data type is not tuple
                if datatype is not tuple:
                    
                    # If the piece of data is available
                    if data[0] != "?":
                        
                        # Convert it to the appropriate data type
                        data = str(datatype(data[0]))
                    
                    # Otherwise
                    else:

                        # Data will be None
                        data = None
                
                # If the data type is tuple   
                else:
                    
                    # Join the pieces of data in a string,
                    # separated by a semicolon
                    data = \
                        ":".join(datatype(\
                            [item for item in data if item != "?"]))

                #-----------------------------------------------------#

                # If the field is the start of the sequence covered by
                # the structure provided in UniProt numbering
                if col == "uniprot_start":

                    # Add the information to the list
                    unp_segment[0] = data

                    # Continue (we do not need to save this field)
                    continue

                #-----------------------------------------------------#

                # If the field is the end of the sequence covered by
                # the structure provided in UniProt numbering
                elif col == "uniprot_end":

                    # Add the information to the list
                    unp_segment[1] = data

                    # Continue (we do not need to save this field)
                    continue

                #-----------------------------------------------------#

                # If the field is the start of the sequence covered by
                # the structure in the structure's own numbering
                elif col == "struct_start":

                    # Add the information to the list
                    struct_segment[0] = data

                    # Continue (we do not need to save this field)
                    continue

                #-----------------------------------------------------#

                # If the field is the end of the sequence covered by
                # the structure in the structure's own numbering
                elif col == "struct_end":

                    # Add the information to the list
                    struct_segment[1] = data

                    # Continue (we do not need to save this field)
                    continue

            #---------------------------------------------------------#
            
            # Otherwise
            else:
                
                # The field is not present
                data = None

            #---------------------------------------------------------#

            # Add the data to the dictionary
            metadata[col] = data

        #-------------------------------------------------------------#

        # Format and add the segment covered by the structure and
        # provided in UniProt numbering to the metadata
        metadata[defaults.DF_COLUMNS["uniprot_segments"]] = \
            defaults.DF_SEPARATORS["range"].join(unp_segment)

        #-------------------------------------------------------------#

        # Format and add the segment covered by the structure and
        # provided in the structure's own numbering to the metadata
        metadata[defaults.DF_COLUMNS["struct_segments"]] = \
            defaults.DF_SEPARATORS["range"].join(struct_segment)

        #-------------------------------------------------------------#

        # Return the dictionary of metadata
        return metadata

