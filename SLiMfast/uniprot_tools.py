#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    uniprot_tools.py
#
#    Utilities to retrieve data from the UniProt database and
#    manipulate them.
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
    "Utilities to retrieve data from the UniProt database and " \
    "manipulate them."


#######################################################################


# Standard library
import logging as log
import re
import requests
# Third-party packages
import pandas as pd


#######################################################################


# Get the module logger
logger = log.getLogger(__name__)


#######################################################################


def read_fasta(fasta_file):
    """Read a FASTA file containing only one sequence and return the
    sequence's header and the sequence itself.

    Parameters
    ----------
    fasta_file : ``str``
        The FASTA file.

    Returns
    -------
    header : ``str``
        The sequence's header.

    sequence : ``str``
        The sequence. 
    """

    # Open the FASTA file
    with open(fasta_file, "r") as f:

        # Initialize the header to an empty string
        header = ""
        
        # Initialize the sequence to an empty string
        sequence = ""

        #-------------------------------------------------------------#
        
        # For each line
        for line in f:
            
            # If you encounter an empty line or a comment line
            if re.match(r"^\s*$", line) or line.startswith("#"):
                
                # Ignore it
                continue

            #---------------------------------------------------------#

            # If the line is the header of the FASTA file
            elif line.startswith(">"):

                # If we have already found a header (= the FASTA file
                # contains multiple sequences)
                if header:

                    # Warn the user
                    warnstr = \
                        f"The FASTA file '{fasta_file}' contains " \
                        "multiple sequences. Only the first " \
                        "one will be read."
                    logger.warning(warnstr)

                    # Break out of the loop
                    break
                
                # Get the header (with the initial '>' stripped)
                header = line.lstrip(">").rstrip("\n")

            #---------------------------------------------------------#

            # If the line is any other line (sequence line)
            else:
                
                # Add the sequence on the line to the full sequence
                # (the sequence may be split on multiple lines)
                sequence += line.rstrip("\n")

        #-------------------------------------------------------------#

        # Return header and sequence
        return header, sequence


def write_uniprot_fasta(uniprot_id,
                        fasta_file_handle):
    """Given a protein's UniProt ID, fetch the protein's sequence
    from UniProt and write it to a FASTA file.

    Parameters
    ----------
    uniprot_id : ``str``
        The protein's UniProt ID.

    fasta_file_handle : ``file handle``
        The FASTA file's handle.
    """

    # The generic URL where the FASTA sequences are found on UniProt
    unp_url = "https://www.uniprot.org/uniprot/{:s}.fasta"

    #-----------------------------------------------------------------#

    # Open the output FASTA file
    with fasta_file_handle as out:
        
        # Get the data from UniProt
        response = requests.get(unp_url.format(uniprot_id))

        # Inform the user that the sequence has been retrieved
        # successfully from UniProt
        infostr = \
            "The FASTA sequence for the protein with ID " \
            f"'{uniprot_id}' was successfully retrieved from " \
            "the UniProt database."
        logger.info(infostr)
        
        # Write the FASTA file
        out.write(response.text)

        # Inform the user that the FASTA file was successfully written
        infostr = \
            "The FASTA sequence for the protein with UniProt ID " \
            f"'{uniprot_id}' was successfully saved to " \
            f"'{fasta_file_handle.name}'."
        logger.info(infostr)


def get_pdbs(uniprot_id):
    """Given a protein's UniProt ID, get the IDs of the structures 
    associated with the protein in the Protein Data Bank.

    Parameters
    ----------
    uniprot_id : ``str``
        The protein's UniProt ID.

    Yields
    ------
    The 
    """

    # The generic URL where the parseable text versions of the data
    # associated to single entries are stores
    unp_url = "https://www.uniprot.org/uniprot/{:s}.txt"

    #-----------------------------------------------------------------#

    # How the line containing the cross-reference with the Protein
    # Data Bank appears in the parseable UniProt file
    CROSS_PDB = "DR   PDB;"

    #-----------------------------------------------------------------#

    # Get the data for the entry associated to the given ID
    response = requests.get(unp_url.format(uniprot_id))

    #-----------------------------------------------------------------#
    
    # For each line in the text file
    for line in response.text.split("\n"):
        
        # If the line corresponds to a cross-reference with the
        # Protein Data Bank 
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

            # Yield the structure's PDB ID
            yield struct_id

    #-----------------------------------------------------------------#

    # Inform the user about the successful retrieval of the structure's
    # IDs
    infostr = \
        "The IDs of structures available in the Protein Data Bank " \
        f"for the protein with UniProt ID '{uniprot_id}' were " \
        "successfully retrieved."
    logger.info(infostr)

