#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    uniprot_tools.py
#
#    Utility functions to retrieve data from the UniProt database.
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
import requests
# Third-party packages
import pandas as pd
# SLiMfast
from SLiMfast.defaults import (
    PDB_CHAIN_FIELD,
    PDB_ID_FIELD,
    PDB_SEGMENTS_FIELD,
    RANGE_SEP,
    SEGMENTS_SEP,
    UNP_ID_FIELD,
    UNP_SEGMENTS_FIELD)
from SLiMfast import pdb_tools



def write_uniprot_fasta(uniprot_id, fasta_file_handle):
    """Fetch the FASTA sequence corresponding to a UniProt
    ID and write it to a FASTA file.
    """

    # Generic URL where the FASTA sequences are found on UniProt
    uniprot_url = "https://www.uniprot.org/uniprot/{:s}.fasta"

    # Open the output FASTA file
    with fasta_file_handle as out:
        # Get the data from the Protein Data Bank
        response = requests.get(uniprot_url.format(uniprot_id))
        # Write the FASTA file
        out.write(response.text)


def get_pdbs(uniprot_id):
    """Get the PDB IDs (and specific chains) associated to
    a UniProt ID.
    """

    # URL where the parseable text versions of the data
    # associated to single entries are stores
    uniprot_url = "https://www.uniprot.org/uniprot/{:s}.txt"

    # Get the data for the entry associated to the given ID
    response = requests.get(uniprot_url.format(uniprot_id))

    # Empty dictionary to store PDB IDs and chains associated
    # to the entry
    pdbs = []
    
    # For each line in the text file
    for line in response.text.split("\n"):
        
        # If the line corresponds to a cross-reference
        # with the Protein Data Bank 
        if line.startswith("DR   PDB;"):
            
            # Split the data into the different fields
            db, pdb_id, exp, res, chain_res = \
                [item.strip(" ") for item \
                 in line.rstrip(".\n").split(";")]
            
            # Get the PDB chain(s) corresponding to the entry
            # pdb_chains, pdb_res = chain_res.split("=")
            # pdb_chains = pdb_chains.split("/")

            # We do not take the mapping between UniProt
            # residues and PDB residues from here because
            # deletions within the chain (so the same PDB
            # segment mapping to multiple UniProt segments)
            # are not reported (i.e. in Q07817)

            # Add the PDB ID, together with all chains that
            # represent the entity and the portion of the
            # chain they cover
            pdbs.append(pdb_id)

    # Return the PDB IDs and chains associated to the UniProt ID
    return pdbs


def get_pdbs_df(uniprot_id):
    """Create a dataframe containing information about the PDB
    structures available for a protein, given the protein's
    UniProt ID.
    """
    
    # Create an empty list where to store the data that
    # will be converted into a dataframe
    data = []

    # Get the PDB IDs corresponding to the given UniProt ID
    pdbs = get_pdbs(uniprot_id = uniprot_id)

    # If no matching PDB structures were found for the
    # UniProt ID
    if not pdbs:
        # Return an empty dataframe
        return pd.DataFrame()
    
    # For each PDB IDcorresponding to the given UniProt ID
    for pdb_id in pdbs:
        
        # Get the mapping between the segments of the protein
        # found in the PDB structure and their corresponding
        # UniProt numbering
        mappings = \
            pdb_tools.get_uniprot_segments(pdb_id = pdb_id,
                                           uniprot_id = uniprot_id)
        
        # For each PDB chain and the corresponding UniProt
        # segments
        for pdb_chain, segments in mappings.items():
            
            # Get the starting and ending points of the UniProt
            # segments
            pdb_segments = \
                [RANGE_SEP.join([str(t) for t in s[0]]) for s in segments]
            unp_segments = \
                [RANGE_SEP.join([str(t) for t in s[1]]) for s in segments]

            pdb_segments_joined = SEGMENTS_SEP.join(pdb_segments)
            unp_segments_joined = SEGMENTS_SEP.join(unp_segments)
            
            # Append the data of interest to the dictionary
            data.append({PDB_ID_FIELD : pdb_id,
                         PDB_CHAIN_FIELD : pdb_chain,
                         UNP_ID_FIELD : uniprot_id,
                         PDB_SEGMENTS_FIELD : pdb_segments_joined,
                         UNP_SEGMENTS_FIELD : unp_segments_joined})

    # Create a dataframe from the data
    df = pd.DataFrame(data)

    # Return the dataframe
    return df
