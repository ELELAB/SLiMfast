#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    util.py
#
#    Utility functions.
#
#    Copyright (C) 2020 Valentina Sora 
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



# standard library
import json
import re
import requests
# third-party packages
import pandas as pd



########################## RUN EXTERNAL TOOLS #########################



def run_hhblits(executable):
    pass


def run_hhsearch(executable):
    pass


def run_dssp(executable):
    pass


def run_naccess(executable):
    pass



################################ OTHERS ###############################



def hhsearch2df(hhsearch_results):
    """Convert data from HHsearch output files to a dataframe.
    """
    

    with open(hhsearch_results, "r") as f:

        # boolean switch to indicate when to parse metadata;
        # on at the beginning because metadata are at the
        # beginning of the file
        parse_metadata = True
        # boolean switch to indicate when to parse the table
        # containing the hits; it comes after the metadata
        parse_hits = False
        # boolean switch to indicate when to parse the alignments;
        # they come after the hits' table
        parse_als = False

        # initialize the indexes of the lines where the hits'
        # table starts and ends
        start_hits_ix = 0
        end_hits_ix = 0

        # initialize the column index where the name of the fields
        # in the metadata end (to separate field names and data)
        metadata_field_end_ix = 14

        # format for data in scientific notation
        sci_fmt = "{:,.1E}"
        
        # set the column-width specifiers (open intervals)
        # for the fields of the hits' table (will be read by
        # pandas as a fixed-width dataframe)
        colspecs_hits = [(0,3), (4,11), (11,34), (35,40), (41,48),
                         (49,57), (58,64), (65,69), (70,74), \
                         (75,84), (85,102)]

        # define the column names for the hits' table
        colnames_hits = ["No", "Hit", "Short description", "Prob", \
                         "E-value", "P-value", "Score", "SS", "Cols", \
                         "Query HMM", "Template HMM"]

        # define the column names for the data retrieved from
        # the aligments
        colnames_als = ["Identities", "Similarity", "Sum_probs", \
                        "Template_Neff", "Organism"]

        # initialize an empty dictionary to store the metadata
        # and one to store the info on the hits gathered from
        # the alignments
        metadata = {}
        info_from_als = {}

        # for each line
        for ix, line in enumerate(f):

            # if you are parsing the metadata section
            if parse_metadata:

                # add the metadata to the dictionary
                metadata[line[0:metadata_field_end_ix].rstrip(" ")] = \
                    line[metadata_field_end_ix].rstrip("\n")

                # if you encounter an empty line
                if re.match(r"^\s*$", line):
                    # turn off metadata parsing,
                    # since the section ended
                    parse_metadata = False
                    # turn on the parsing for the hits' table
                    # since it is the next section
                    parse_hits = True
                    # the starting index for the hits' table
                    # section is the one of the next line
                    start_hits_ix = ix + 1
                    # ignore the current (empty) line
                    continue

            # if you are parsing the hits' table section
            if parse_hits:

                # if you encounter a blank line
                if re.match(r"^\s*$", line):
                    # turn off the parsing of the hits' table
                    parse_hits = False
                    # turn on the parsing for the alignments
                    parse_als = True
                    # the ending index for the hits' table
                    # is the one of the current line
                    # (the interval is open, so the table
                    # ends on the previous line)
                    end_hits_ix = ix
                    # ignore the current (empty) line
                    continue

            # if you are parsing the alignments' section
            if parse_als:
                
                # if the line is the header line 
                if line.startswith("No "):
                    # get the hit number
                    hit_num = line.rstrip("\n").lstrip("No ")
                
                # if the line is the first line of the metadata
                # associated with the alignment
                elif line.startswith(">"):
                    # find the name of the hit 
                    hit_name = line.split(" ")[0].lstrip(">")
                    hit_name = "{0}.{1}".format(\
                                hit_name.split("_")[0].lower(), \
                                hit_name.split("_")[1])
                    
                    # find the name of the organism with a regular
                    # expression and get rid of the parentheses
                    org = re.findall(r"{.*?}", \
                            line)[0].rstrip("}").lstrip("{").lower()
                    org = org[0].upper() + org[1:]

                # if the line is the second line of the metadata
                # associated with the alignment
                elif line.startswith("Probab"):
                    # get all the data on it
                    prob, evalue, score, aligned_cols, \
                    identities, similarity, sum_probs, \
                    template_neff = \
                        [float(item.split("=")[1].rstrip("%")) \
                         for item in line.split(" ") \
                         if item != ""]
                    # data to add to the corresponding hit in the table
                    data_to_add = [identities, similarity, sum_probs, \
                                   template_neff, org]
                    # generate a hit ID compatible with the one used
                    # in the table
                    hit_id = f"{hit_num}_{hit_name}"
                    # add the info from the alignment to the dictionary
                    # storing them for all the hits
                    info_from_als[hit_id] = \
                        dict(colnames_als, data_to_add)

        # read the 'hit' section as a fixed-width-formatted file
        # into a dataframe
        hits_df = pd.read_fwf(filename, \
                              skiprows = start_hits_ix, \
                              nrows = end_hits_ix-start_hits_ix-1, \
                              header = 0, \
                              colspecs = colspecs_hits)

        # specify row and column names
        hits_df.columns = colnames_hits
        
        # row names are composed by both the hit number and the
        # PDB ID (with chain) to both ensure their uniqueness and
        # be informative at the same time
        hits_df.index = [str(num) + "_" + hit for num, hit \
                         in zip(hits_df["No"], hits_df["Hit"])]

        # add two columns storing just the PDB ID and the chain
        hits_df[["PDB_Id", "Chain"]] = hits_df["Hit"].str.split("_", 1)

        # convert information from the alignments to a dataframe
        als_df = pd.DataFrame.from_dict(info_from_als, \
                                        orient = "index")

        # add the information to the original dataframe
        complete_df = pd.concat([hits_df, als_df], \
                                 axis = 1, \
                                 sort = False)

        # specify the format for data in scientific notation
        complete_df["E-value"] = \
            complete_df["E-value"].map(sci_fmt.format)
        complete_df["P-value"] = \
            complete_df["P-value"].map(sci_fmt.format)

        # return the dataframe
        return complete_df


def pdbchain2uniprot(pdb_id, pdb_chain):
    """Get the UniProt ID associated to a chain in a PDB.
    """

    # URL where the mapping between PDB entities in a PDB entry
    # and their corresponding UniProt data is stored on PDBe
    mapping_url = \
        "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"

    # name of the field storing UniProt data
    up_field = "UniProt"
    # name of the field storing the chain ID
    chain_id_field = "chain_id"
    # name of the field storing the mappings
    mappings_field = "mappings"

    # get the data from PDBe
    response = requests.get(mapping_url.format(pdb_id))
    # convert to text and read the resulting dictionary through JSON
    data = json.loads(response.text)

    # convert the PDB ID in all lowercase since it it stored
    # like this in the mapping file
    pdb_id_low = pdb_id.lower()

    # for each UniProt ID associated to an entity in the PDB
    # and associated data
    for up_id, up_data in data[pdb_id_low][up_field].items():
        # if the entity with such UniProt ID is the one with
        # the desided chain ID
        if up_data[mappings_field][0][chain_id_field] == chain_id:
            # return the UniProt ID
            return up_id

    # raise an error if not UniProt ID was found for the
    # desired chain
    errmsg = f"UniProt ID not found for PDB ID " \
             f"{pdb_id} chain {chain_id}."
    
    raise ValueError(errmsg)


