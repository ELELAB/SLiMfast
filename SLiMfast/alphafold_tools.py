#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    alphafold_tools.py
#
#    Utility functions to retrieve data from PDB/mmCIF files from
#    the AlphaFold Protein Structure Database.
#
#    Copyright (C) 2022 Valentina Sora 
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
import itertools
import json
import logging as log
import requests
from tempfile import NamedTemporaryFile
# Third-party packages
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
# SLiMfast
from SLiMfast.defaults import (
    ALPHAFOLD_ID_FIELD,
    DEPOSITION_DATE_FIELD,
    RANGE_SEP,
    STRUCT_CHAIN_FIELD,
    STRUCT_ID_FIELD,
    STRUCT_SEGMENTS_FIELD,
    UNP_ISOFORM_FIELD,
    UNP_SEGMENTS_FIELD,
    )



# Get the module logger
logger = log.getLogger(__name__)



# Metadata (taken from MMCIF files) retrieved for AlphaFold
# IDs of interest
STRUCTURE_METADATA = [
    
    # Beginning of the sequence covered by the structure
    # in the UniProt numbering
    ("_ma_target_ref_db_details.seq_db_align_begin",
     int,
     "uniprot_start"),
    # End of the sequence covered by the structure
    # in the UniProt numbering
    ("_ma_target_ref_db_details.seq_db_align_end",
     int,
     "uniprot_end"),
    # Begin of the sequence covered by the structure
    # in the PDB numbering
    ("_struct_ref.pdbx_align_begin",
     int,
     "struct_start"),
    # End of the sequence covered by the structure
    # in the PDB numbering
    ("_struct_ref.pdbx_align_begin",
     int,
     "struct_end"),
    # AlphaFold ID
    ("_database_2.database_code",
     str, 
     ALPHAFOLD_ID_FIELD),
    # Chain ID
    ("_entity_poly.pdbx_strand_id",
     str,
     STRUCT_CHAIN_FIELD),
    # UniProt isoform
    ("_ma_target_ref_db_details.seq_db_align_isoform",
     str,
     UNP_ISOFORM_FIELD),
    # Deposition date
    ("_pdbx_database_status.recvd_initial_deposition_date",
     str,
     DEPOSITION_DATE_FIELD),
    ]



def get_alphafold_data(uniprot_id = None,
                       struct_id = None):
    """Get the PDB structure(s) associated to a specific
    UniProt ID in the AlphaFold Protein Structure Database,
    or a specific structure from the database.
    """

    # Check if either the UniProt ID or the structure ID have
    # been passed
    if uniprot_id is None and struct_id is None:

        # Raise an error
        errstr = \
            "You must provide either the UniProt ID (to find " \
            "all structures associated with it in the AlphaFold " \
            "Protein Structure Database) or the specific structure " \
            "ID (i.e., AF-P10415-F1-model_v3)."
        raise ValueError(errstr)

    # Generic string representing the ID of the
    # structures deposited in AlphaFoldDB
    struct_id_fmt = "AF-{:s}-F{:d}-model_v{:d}"

    # If the structure ID has been specified
    if struct_id is not None:

        # This will become the string to be used to
        # find the structure in the database
        struct_id_fmt = struct_id 
    
    # Generic URL where AlphaFoldDB PDB files are found
    struct_url = \
        f"https://www.alphafold.ebi.ac.uk/files/" \
        f"{struct_id_fmt}.pdb"


    #-------------------- Structure ID specified ---------------------#


    # If the structure ID was specified
    if struct_id is not None:

        # Try to get the data from the AlphaFold Protein Structure
        # Database
        response = requests.get(struct_url)

        # If something went wrong
        if not response.ok:

            # Re-raise the exception
            response.raise_for_status()

        # Return the data
        return [l for l in response.text.split("\n")]


    #------------------ Structure ID not specified -------------------#


    # Create an empty dictionary to store the structures
    alphafold_data = {}

    # Try to get all versions of the structure
    for version in itertools.count(start = 1):

        # Try to get all fragments (if the structure
        # is split in multiple fragments
        for fragment in itertools.count(start = 1):
            
            # Try to get the data from the AlphaFold
            # Protein Structure Database
            response = \
                requests.get(struct_url.format(uniprot_id,
                                               fragment,
                                               version))

            # If the structure was not found (it does
            # not exist)
            if response.status_code == 404:

                # Break out of the loop since we
                # do not have more fragments
                break

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

        # If the structure was not found and it is because
        # we have reached the latest version
        if response.status_code == 404 and fragment == 1:

            # Break out of the loop since we
            # do not have more versions
            break

    # Return the data retrieved
    return alphafold_data


def get_structure_metadata(struct_id):
    """Get metadata associated to a structure in the 
    AlphaFoldDB from its mmCIF file.
    """

    # Metadata are collected from mmCIF files and not
    # directly from the PDB files since the mmCIF format
    # is more consistent and contains more information

    # Generic URL where mmCIF files for single entries
    # are found
    mmcif_url = \
        "https://www.alphafold.ebi.ac.uk/files/{:s}.cif"
    
    # Create a dictionary to store the metadata
    metadata = {}

    # Use a temporary file to store the mmCIF file
    with NamedTemporaryFile("w") as tmp:
        
        # Get the mmCIF file
        response = requests.get(mmcif_url.format(struct_id))

        # Write the mmCIF file to the temporary file
        tmp.write(response.text)
         
        # Create a dictionary from the header data
        mmcif_dict = MMCIF2Dict(tmp.name)

        # Initialize a list to store the segment
        # covered by the AlphaFold structure in the
        # UniProt numbering
        unp_segment = [None, None]

        # Initialize a list to store the segment
        # covered by the AlphaFold structure in the PDB
        # numbering
        struct_segment = [None, None]

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

                        # Data will be None
                        data = None
                
                # If the data type is tuple   
                else:
                    
                    # Join the pieces of data in a string,
                    # separated by a semicolon
                    data = \
                        ":".join(datatype(\
                            [item for item in data if item != "?"]))


                # If the field is the start of the sequence
                # covered by the structure in the UniProt
                # numbering
                if col == "uniprot_start":

                    # Add the information to the list
                    unp_segment[0] = data

                    # Continue (we do not need to save
                    # this field)
                    continue

                # If the field is the end of the sequence
                # covered by the structure in the UniProt
                # numbering
                elif col == "uniprot_end":

                    # Add the information to the list
                    unp_segment[1] = data

                    # Continue (we do not need to save
                    # this field)
                    continue

                # If the field is the start of the sequence
                # covered by the structure in the PDB
                # numbering
                elif col == "struct_start":

                    # Add the information to the list
                    struct_segment[0] = data

                    # Continue (we do not need to save
                    # this field)
                    continue

                # If the field is the end of the sequence
                # covered by the structure in the PDB
                # numbering
                elif col == "struct_end":

                    # Add the information to the list
                    struct_segment[1] = data

                    # Continue (we do not need to save
                    # this field)
                    continue
            
            # Otherwise
            else:
                
                # The field is not present
                data = None

            # Add data to the dictionary
            metadata[col] = data

        # Format and add the sequence segment covered
        # by the structure in the UniProt numbering
        # to the metadata
        metadata[UNP_SEGMENTS_FIELD] = \
            RANGE_SEP.join(unp_segment)

        # Format and add the sequence segment covered
        # by the structure in the PDB numbering
        # to the metadata
        metadata[STRUCT_SEGMENTS_FIELD] = \
            RANGE_SEP.join(struct_segment)

        # Return the dictionary of metadata
        return metadata