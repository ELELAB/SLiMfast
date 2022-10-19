#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    test.py
#
#    Unit tests.
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
from collections import defaultdict
import filecmp
import json
import os
import os.path as op
from tempfile import NamedTemporaryFile
# Third-party packages
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
# SLiMfast
from SLiMfast import (
    disprot_tools,
    dssp_tools,
    hhsearch_tools,
    pdb_tools,
    uniprot_tools,
    util)



########################## HELPER FUNCTIONS ###########################



def load_dssp_csv(csv_file):

    # Get the data types to be used for all colums of the dataframe
    dtype = dssp_tools.DSSPParser.DTYPES

    # Load the data frame
    return pd.read_csv(csv_file, index_col = 0, dtype = dtype)


def load_hhsearch_csv(csv_file):

    # Get the data types to be used for all colums of the dataframe
    dtype = hhsearch_tools.HHsearchParser.DTYPES

    # Load the data frame
    return pd.read_csv(csv_file, index_col = 0, dtype = dtype)


def load_json(json_file):

    # Load a JSON file content as a dictionary
    return json.load(open(json_file, "r"))


def load_data(txt_file):
    
    # Load lines from a text file
    return [l.rstrip("\n") for l in open(txt_file, "r").readlines()]


def list_dir(top_dir):

    # List all top level non-hidden directories in a given directory
    return [d for d in next(os.walk(op.join(top_dir)))[1] \
            if not d.startswith(".")]



############################### FIXTURES ##############################



@pytest.fixture
def data_dir(request):

    # Directory where data for tests are stored
    return op.join(request.fspath.dirname, "data")


@pytest.fixture
def examples(data_dir):
    
    # Create an empty defaultidict of defaultdicts to store data
    # about the examples
    data = []

    # For each directory containing an example (which has a
    # UniProt ID as name)
    for uniprot_id in list_dir(data_dir):

        e_data = defaultdict(lambda : defaultdict(dict))
        
        # Get the path to the current example's directory
        e_dir = op.join(data_dir, str(uniprot_id))

        # Load data regarding HHsearch-related tests
        e_data["hhsearch"] = \
          {"test" : op.join(e_dir, "hhsearch.hhr"),
           "df" : load_hhsearch_csv(op.join(e_dir, "hhsearch.csv")),
           "als" : load_json(op.join(e_dir, "hhsearch.json")),
           "filtered_df" : \
             load_hhsearch_csv(op.join(e_dir, "hhsearch_filt.csv")),
           "sorted_df" : \
             load_hhsearch_csv(op.join(e_dir, "hhsearch_sort.csv"))}

        # Get the UniProt ID and the FASTA file containing
        # the protein sequence
        e_data["uniprot"] = \
            {"id" : uniprot_id, 
             "fasta" : op.join(e_dir, "sequence.fasta")}

        # Load data regarding PDB-related tests
        e_data["pdb"] = \
          [{"id" : p, 
            "segments" : \
              load_json(op.join(e_dir, f"{p}/segments.json")),
            "data" : \
              load_data(op.join(e_dir, f"{p}/structure.pdb")),
            "miss_res" : \
              load_json(op.join(e_dir, f"{p}/miss_res.json")),
            "metadata" : \
              load_json(op.join(e_dir, f"{p}/metadata.json")),
            "dssp" : \
              {"test" : op.join(e_dir, f"{p}/dssp.out"),
               "df" : load_dssp_csv(op.join(e_dir, f"{p}/dssp.csv")),
               "seq" : \
                 load_json(op.join(e_dir, f"{p}/dssp_seq.json")),
               "seqdssp" : \
                 load_json(op.join(e_dir, f"{p}/dssp_seqdssp.json"))}}
             for p in list_dir(e_dir)]

        # Save the data for the single example        
        data.append(e_data)

    # Return data for all examples
    return data



################################ TESTS ################################



def test_disprot_parser(examples):
    
    # Create the DisProt annotations parser
    parser = disprot_tools.DisProtParser()
    
    # For each example
    for example in examples:

        # Test the creation of the dataframe of DisProt annotations
        df = parser.get_annotations(example["uniprot"]["id"])


def test_dssp_parser(examples):

    # Create the DSSP output parser
    parser = dssp_tools.DSSPParser()
    
    # For each example
    for example in examples:

        # For each PDB
        for pdb in example["pdb"]:

            # Get the dssp data
            data = pdb["dssp"]
        
            # Get the dataframe, the FASTA sequence and the
            # FASTA-like DSSP sequence from the sample data
            df, seq, seqdssp = parser.parse(data["test"])

            # Check that the results are equal to the expected ones
            assert_frame_equal(df, data["df"])
            assert seq == data["seq"]
            assert seqdssp == data["seqdssp"]


def test_hhsearch_parser(examples):

    # Set the filters and sorting parameters to be used for
    # the output dataframes
    filters = {"min_uniprot_coverage_percent" : 70,
               "min_template_neff" : 5.0}
    sort_by = ["uniprot_coverage_percent", "similarity"]
    ascending = [False, False]

    # Create the HHsearch output parser
    parser = hhsearch_tools.HHsearchParser()
    
    # For each example
    for example in examples:
        
        # Get the dataframe containing the hits and the
        # dictionary containing the alignments 
        df, als = parser.parse(example["hhsearch"]["test"])
        
        # Filter the dataframe
        filt_df = hhsearch_tools.filter_hits(df, filters)
        
        # Sort the dataframe
        sort_df = hhsearch_tools.sort_hits(df, sort_by, ascending)

        # Check that the results are equal to the expected ones
        assert_frame_equal(df, example["hhsearch"]["df"])
        assert_frame_equal(filt_df, example["hhsearch"]["filtered_df"])
        assert_frame_equal(sort_df, example["hhsearch"]["sorted_df"])
        assert als == example["hhsearch"]["als"]


def test_get_uniprot_segments(examples):

    # For each example  
    for example in examples:

        # Get the UniProt ID
        uniprot_id = example["uniprot"]["id"]

        # For each PDB
        for pdb in example["pdb"]:
            
            # Get the PDB ID
            pdb_id = pdb["id"]

            # Get the protein segments covered by the PDB entry with
            # that ID
            segs = pdb_tools.get_protein_segments(pdb_id, uniprot_id)

            # Check that the results are equal to the expected ones
            assert segs == pdb["segments"]


def test_get_pdb_data(examples):
    
    # For each example
    for example in examples:

        # For each PDB
        for pdb in example["pdb"]:

            # Get the PDB ID
            pdb_id = pdb["id"]

            # Get the content of the PDB entry having that ID
            data = pdb_tools.get_pdb_data(pdb_id)

            # Check that the results are equal to the expected ones
            assert data == pdb["data"]


def test_get_missing_residues(examples):

    # For each example
    for example in examples:

        # For each PDB
        for pdb in example["pdb"]:

            # Get the content of the PDB entry
            pdb_data = pdb["data"]

            # Get the missing residues from the content of the
            # entry
            miss_res = pdb_tools.get_missing_residues(pdb_data)

            # Check that the results are equal to the expected ones
            assert miss_res == pdb["miss_res"]


def test_get_structure_metadata(examples):

    # For each example
    for example in examples:

        # For each PDB
        for pdb in example["pdb"]:

            # Get the PDB ID
            pdb_id = pdb["id"]

            # Get the metadata of the structure with that ID
            metadata = pdb_tools.get_structure_metadata(pdb_id)

            # Check that the results are equal to the expected ones
            assert metadata == pdb["metadata"]

