#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    uniprot_tools.py
#
#    Utility functions to retrieve data from the UniProt database.
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
import logging as log
import requests
# Third-party packages
import pandas as pd



# Get the module logger
logger = log.getLogger(__name__)



def write_uniprot_fasta(uniprot_id,
                        fasta_file_handle):
    """Fetch the FASTA sequence corresponding to a UniProt
    ID and write it to a FASTA file.
    """

    # Generic URL where the FASTA sequences are found on UniProt
    UNP_URL = "https://www.uniprot.org/uniprot/{:s}.fasta"

    # Open the output FASTA file
    with fasta_file_handle as out:
        
        # Get the data from UniProt
        response = requests.get(UNP_URL.format(uniprot_id))

        # Inform the user that the sequence has been retrieved
        # successfully from UniProt
        infostr = \
            f"The FASTA sequence for the protein with ID " \
            f"{uniprot_id} has been successfully retrieved from " \
            f"the UniProt database."
        logger.info(infostr)
        
        # Write the FASTA file
        out.write(response.text)

    # Inform the user that the FASTA file has been written
    # successfully
    infostr = \
        f"The FASTA sequence for the protein with UniProt ID " \
        f"{uniprot_id} has been successfully written to " \
        f"{fasta_file_handle.name}."
    logger.info(infostr)


def get_pdbs(uniprot_id):
    """Get the IDs (and specific chains) associated to
    a UniProt ID in the Protein Data Bank.
    """

    # URL where the parseable text versions of the data
    # associated to single entries are stores
    UNP_URL = "https://www.uniprot.org/uniprot/{:s}.txt"

    # How the line containing the cross-reference with the
    # Protein Data bank appears in the UniProt entry
    CROSS_PDB = "DR   PDB;"

    # Get the data for the entry associated to the given ID
    response = requests.get(UNP_URL.format(uniprot_id))

    # Empty dictionary to store PDB IDs and chains associated
    # to the entry
    pdbs = []
    
    # For each line in the text file
    for line in response.text.split("\n"):
        
        # If the line corresponds to a cross-reference
        # with the Protein Data Bank 
        if line.startswith(CROSS_PDB):
            
            # Split the data into the different fields
            db, struct_id, exp, res, chain_res = \
                [item.strip(" ") for item \
                 in line.rstrip(".\n").split(";")]

            # We do not take the mapping between UniProt
            # residues and PDB residues from here because
            # deletions within the chain (so the same PDB
            # segment mapping to multiple UniProt segments)
            # are not reported (i.e. in Q07817)

            # Add the PDB ID, together with all chains that
            # represent the entity and the portion of the
            # chain they cover
            pdbs.append(struct_id)

    # Inform the user about the successful creation of the
    # data frame with the available PDB structures
    infostr = \
        f"The data frame containing the available structures " \
        f"for the protein with UniProt ID {uniprot_id} has been " \
        f"successfully created."
    logger.info(infostr)

    # Return the PDB IDs and chains associated to the UniProt ID
    return pdbs
