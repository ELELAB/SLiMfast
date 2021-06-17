#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    util.py
#
#    Miscellanea utility functions for the SLiMfast pipeline.
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
import argparse
import logging as log
import math
import os
import re
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
    STRUCTURE_METADATA,
    UNP_ID_FIELD,
    UNP_SEGMENTS_FIELD)
from SLiMfast import (
    disprot_tools,
    dssp_tools,
    hhsearch_tools,
    pdb_tools,
    uniprot_tools)



# Get the module logger
logger = log.getLogger(__name__)


############################## CONFIGURE ##############################



def get_arguments(description,
                  log_file_default = None,
                  **kwargs):
    """Create the argument parser, add the arguments requested
    and return them.
    """

    # Create the argument parser
    parser = argparse.ArgumentParser(description = description)
    
    # Add arguments
    if "uniprot_id" in kwargs.keys():
        id_help = "UniProt ID."
        parser.add_argument("-id", "--uniprot-id",
                             type = str,
                             required = True,
                             help = id_help)

    if "input_protein_file" in kwargs.keys():
        ip_help = \
            "Input CSV file containing the available PDB " \
            "structures of the protein."
        parser.add_argument("-ip", "--input-protein-file",
                             type = str,
                             required = True,
                             help = ip_help)

    if "input_homologs_file" in kwargs.keys():
        ih_help = \
            "Input CSV file containing the available PDB " \
            "structures of the homologs."
        parser.add_argument("-ih", "--input-homologs-file",
                             type = str,
                             required = True,
                             help = ih_help)

    if "input_alignments_file" in kwargs.keys():
        ia_help = "Name of the input JSON alignments file."
        parser.add_argument("-ia", "--input-alignments-file",
                             type = str,
                             required = True,
                             help = ia_help)

    if "input_slims_file" in kwargs.keys():
        is_help = "SLiMs in the UniProt sequence."
        parser.add_argument("-is", "--input-slims-file",
                            type = str,
                            required = True,
                            help = is_help)

    if "output_protein_file" in kwargs.keys():
        op_default = "protein.csv"
        op_help = \
            f"Name of the output CSV file. Default is " \
            f"'{op_default}'."
        parser.add_argument("-op", "--output-protein-file",
                             type = str,
                             default = op_default,
                             help = op_help)

    if "output_homologs_file" in kwargs.keys():
        oh_default = "homologs.csv"
        oh_help = \
            f"Name of the output CSV file. Default is " \
            f"'{oh_default}'."
        parser.add_argument("-oh", "--output-homologs-file",
                             type = str,
                             default = oh_default,
                             help = oh_help)

    if "output_alignments_file" in kwargs.keys():
        oa_default = "alignments.json"
        oa_help = \
            f"Name of the output JSON alignments file. Default is " \
            f"'{oa_default}'."
        parser.add_argument("-oa", "--output-alignments-file",
                             type = str,
                             default = oa_default,
                             help = oa_help)

    if "output_slims_protein_file" in kwargs.keys():
        osp_default = "slims_protein.csv"
        osp_help = \
            f"Name of the output CSV file containing SLiMs found " \
            f"in available structures of the protein. Default is " \
            f"'{osp_default}'."
        parser.add_argument("-osp", "--output-slims-protein-file",
                             type = str,
                             default = osp_default,
                             help = osp_help)

    if "output_slims_homologs_file" in kwargs.keys():
        osh_default = "slims_homologs.csv"
        osh_help = \
            f"Name of the output CSV file containing SLiMs found " \
            f"in available structures of the homologs. Default is " \
            f"'{osh_default}'."
        parser.add_argument("-osh", "--output-slims-homologs-file",
                             type = str,
                             default = osh_default,
                             help = osh_help)

    if "output_slims_disprot_file" in kwargs.keys():
        osd_default = "slims_disprot.csv"
        osd_help = \
            f"Name of the output CSV file containing SLiMs not found " \
            f"in the available structures either of the protein or " \
            f"of the homologs, but annotated as disordered in " \
            f"DisProt. Default is '{osd_default}'."
        parser.add_argument("-osd", "--output-slims-disprot-file",
                             type = str,
                             default = osd_default,
                             help = osd_help)

    if "configrun_file" in kwargs.keys():
        cr_help = "YAML configuration file for running."
        parser.add_argument("-cr", "--configrun-file",
                            type = str,
                            required = True,
                            help = cr_help)

    if "log_file" in kwargs.keys():
        lf_default = log_file_default
        lf_help = f"Name of the log file. Default is {lf_default}."
        parser.add_argument("-lf", "--log-file",
                            type = str,
                            default = lf_default,
                            help = lf_help)
    if "v" in kwargs.keys():
        v_help = "Verbose logging."
        parser.add_argument("-v", "--verbose",
                            action = "store_true",
                            help = v_help)

    if "vv" in kwargs.keys():
        vv_help = "More verbose (DEBUG) logging."
        parser.add_argument("-vv", "--verbose-debug",
                            action = "store_true",
                            help = vv_help)

    # Parse the command line arguments
    args = parser.parse_args()

    # Return the arguments
    return args



############################## READ FILES #############################



def read_fasta(fasta_file):
    """Read header and sequence from a FASTA file.
    """

    with open(fasta_file, "r") as f:

        # Initialize the header to an empty string
        header = ""
        # Initialize the sequence to an empty string
        sequence = ""
        
        # For each lune
        for line in f:
            
            # If you encounter an empty line
            if re.match(r"^\s*$", line):
                # Ignore it
                continue

            # If the line is the header of the FASTA file
            elif line.startswith(">"):
                # Get the header
                header = line.lstrip(">").rstrip("\n")

            # If the line is any other line (sequence line)
            else:
                # Add the sequence on the line to the full sequence
                # (the sequence may be split on multiple lines)
                sequence += line.rstrip("\n")

        # Return header and sequence
        return header, sequence


def read_slims_file(slims_file):
    """Get the instances of a SLiM in a sequence
    from a text file.
    """

    with open(slims_file, "r") as f:

        # Empty list to store SLiM instances
        slim_instances = []

        for line in f:
            # Ignore empty lines and comment lines
            if re.match(r"^\s*$", line) or line.startswith("#"):
                continue
            else:
                slim_instances.append(\
                    tuple(line.rstrip("\n").split(",")))

        # Return the list of SLiM instances
        return slim_instances



############################ GET STRUCTURES ###########################



def add_structures_metadata(df):
    """Add metadata regarding each structure to the
    dataframe.
    """

    # Create an empty list to store the metadata
    # for each unique PDB ID
    metadata_list = []

    # For each unique PDB ID
    for pdb_id in df[PDB_ID_FIELD].unique():
        
        # Try to get the structure metadata for the PDB
        try:
            metadata = \
                pdb_tools.get_structure_metadata(pdb_id = pdb_id)
        
        # If something went wrong, log it and skip the structure
        except Exception as e:
            warnstr = \
                f"Could not fetch metadata for the PDB structure " \
                f"with ID {pdb_id}. This structure will be " \
                f"ignored. Exception: {e}"
            # Log a warning
            logger.warning(warnstr)
            # Go to the next structure
            continue

        # Add a field for the PDB ID (it is needed to merge
        # the dataframes later)
        metadata.update({PDB_ID_FIELD : pdb_id})
        # Add the dictionary of metadata 
        metadata_list.append(metadata)

        # Inform the user that the metadata have been collected
        # successfully
        infostr = f"Metadata collected for PDB {pdb_id}"
        logger.info(infostr)

    # Create a dataframe containing the metadata for each PDB
    # structure
    metadata_df = pd.DataFrame(metadata_list)

    # Combine the concatenated dataframe with the dataframe
    # containing the structures' metadata
    updated_df = pd.merge(left = df,
                          right = metadata_df,
                          on = PDB_ID_FIELD,
                          how = "outer")

    # Return the updated dataframe
    return updated_df


def get_structures_data_protein(uniprot_id,
                                **kwargs):
    """Get the data regarding the available PDB structures
    of the protein.
    """


    #---------------- Get the available PDB structures ---------------#


    # Try to get the dataframe of available PDB structures
    try:
        df = uniprot_tools.get_pdbs_df(uniprot_id = uniprot_id)
    
    # If something went wrong, log it and raise an exception
    except Exception as e:
        errstr = \
            f"Could not get the available PDB structures for the " \
            f"protein with UniProt ID {uniprot_id}. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)


    #------------------------ Add the metadata -----------------------#


    # Try to add metadata for the structures
    try:
        df = add_structures_metadata(df = df)
    
    # If something went wrong, log it and raise an exception
    except Exception as e:
        errstr = \
            f"Could add metadata to the available PDB structures " \
            f"for the protein with UniProt ID {uniprot_id}. " \
            f"Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)


    #----------------------- Sort the datarame -----------------------#


    # Try to sort the dataframe according to the given criteria
    try:
        df = df.sort_values(by = kwargs["sort"]["sort_by"], 
                            ascending = kwargs["sort"]["ascending"])
 
     # If something went wrong, log it and raise an exception   
    except Exception as e:
        errstr = \
            f"Could not sort the dataframe of the available PDB " \
            f"structures for the protein with UniProt ID " \
            f"{uniprot_id}. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)

    # Inform the user that all data has been collected
    # successfully
    infostr = \
        f"All PDB data for {uniprot_id} have been collected " \
        f"successfully."
    logger.info(infostr)

    # Return the dataframe
    return df


def get_structures_data_homologs(uniprot_id,
                                 hhblits_params,
                                 hhsearch_params,
                                 wd,
                                 **kwargs):
    """Get the data regarding the PDB structures of the
    homologs.
    """


    #--------------------------- Run HHpred --------------------------#


    # Try to run HHpred
    try:
        _, hhsearch_output = \
            run_tools.run_hhpred(uniprot_id = uniprot_id,
                                 hhblits_params = hhblits_params,
                                 hhsearch_params = hhsearch_params,
                                 wd = wd)

    # If something went wrong, log it and raise an exception
    except Exception as e:
        errstr = f"Could not complete the HHpred run. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)


    #------------------ Get the homologs' structures -----------------#


    # Create the parser to parse the HHsearch output
    parser = hhsearch_tools.HHsearchParser()

    # Try to parse the HHsearch output to get the list of homologs
    # of interest and their associated data
    try:
        hits, als = parser.parse(hhsearch_output = hhsearch_output)

    # If something went wrong, log it and raise an exception
    except Exception as e:
        errstr = \
            f"Could not get the available PDB structures for the " \
            f"homologs of the protein with UniProt ID " \
            f"{uniprot_id}. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)


    #------------------------ Add the metadata -----------------------#


    # Try to add structures' metadata
    try:
        hits = add_structures_metadata(df = hits)
    
    # If something went wrong, log it and raise an exception
    except Exception as e:
        errstr = \
            f"Could add metadata to the available PDB structures " \
            f"of the homologs of the protein with UniProt ID " \
            f"{uniprot_id}. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)


    #---------------------- Filter the dataframe ---------------------#


    # Try to filter the hits
    try:
        hits = hhsearch_tools.filter_hits(hits = hits,
                                          **kwargs["filters"])
    
    # If something went wrong, log it and raise an exception
    except Exception as e:
        errstr = \
            f"Could not filter the dataframe of the available PDB " \
            f"structures for the homologs of the protein with " \
            f"UniProt ID {uniprot_id}. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)


    #----------------------- Sort the datarame -----------------------#


    # Try to sort the hits
    try:
        hits = hhsearch_tools.sort_hits(hits = hits, 
                                        **kwargs["sort"])
    
    # If something went wrong, log it and raise an exception
    except Exception as e:
        errstr = \
            f"Could not sort the dataframe of the available PDB " \
            f"structures for the homologs of the protein with " \
            f"UniProt ID {uniprot_id}. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)

    # Return the datafame containing the data about the homologs
    # and the dictionary of alignments
    return hits, als



############################## FIND SLIMS #############################



def region_in_protein_segments(seq_range,
                               uniprot_segments):
    """Check if a sequence belongs at least to one UniProt
    segment, given the sequence numeric range.
    """
    
    # Get starting and ending point of the sequence range  
    start_s, end_s = seq_range

    # Get the length of the sequence
    length = len(range(start, end))

    # Get the numeric ranges of all UniProt segments passed
    uniprot_ranges = \
        [(int(s.split(RANGE_SEP)[0]), int(s.split(RANGE_SEP)[1])) \
         for s in uniprot_segments.split(SEGMENTS_SEP)]

    # For each UniProt segment
    for uniprot_start, uniprot_end in uniprot_ranges:

        # If the SLiM is found within the segment
        if start_s >= uniprot_start and end_s <= uniprot_end:
            # Return True                   
            return True

    # If the sequence was never found, return False
    return False


def region_in_homolog_segments(pdb_data,
                               pdb_id,
                               pdb_chain,
                               hit_al,
                               seq_range):
    """Check if a region of the protein of interest is found in the
    PDB segments corresponding to the sequences aligned to such a
    region in the proteins' homolog.
    """
    

    #-------------------- Find the aligned sequence ------------------#


    # Try to find the sequence in the homolog aligned to the
    # sequence of interest in the HHsearch alignments
    try:
        start_in_t, slim_in_t, end_in_t, is_match = \
            hhsearch_tools.get_aligned_sequence(hit_al = hit_al,
                                                seq_range = seq_range)
    
    # If something went wrong, log it and raise an exception
    except Exception as e:
        start, end = seq_range
        errstr = \
            f"Could not find the sequence {start}-{end} in the " \
            f"alignment. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)

    # If the aligned sequence does not perfectly match the sequence
    # of interest, return
    if not is_match:
        return False, None, None, None

    # Get the numeric range covered by the aligned sequence in
    # the PDB structure of the homolog
    slim_range = set(list(range(start_in_t, end_in_t)))


    #------------------- Find the missing residues -------------------#


    # Try to get the missing residues in the PDB structure
    try:
        missing_res = \
            pdb_tools.get_missing_residues(pdb_data = pdb_data)
    
    # If something went wrong, log it and raise an exception
    except Exception as e:
        errstr = \
            f"Could not get the missing residues for the PDB " \
            f"structure with ID {pdb_id}. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)
    
    # Get the missing residues for the chain of interest
    missing_res_chain = set(missing_res[pdb_chain])

    # If none of the residues in the aligned sequence of interest
    # is missing in the structure of the homolog
    if not set(slim_range).intersection(missing_res_chain):
        return True, start_in_t, slim_in_t, end_in_t

    # Return None for each of the structure-related info if
    # some residues of the aligned sequence of interest were
    # missing in the structure of the homolog
    return True, None, None, None


def find_slims(source,
               structures_df,
               slims_to_find,
               find_in_pdb,
               wd,
               als = None):
    """Find the SLiMs in the structures of the protein
    or in those of the homologs.
    """
        
    # Create an empty dictionary storing the data
    # for the PDB structures of interest
    pdb_data_dict = {}

    # Create an empty list to store the PDB files
    # and the SLiMs to be checked for disorder
    structures_and_slims = []

    # Create a copy of the list of SLiMs to be found
    # (so that the original one will not be modified)
    still_to_find = copy.deepcopy(slims_to_find)

    for index, row in structures_df.iterrows():

        # Create an empty list to store the SLiMs
        # found in each PDB
        slims_in_pdb = []

        # Get PDB ID, chain, the segments of the protein
        # of interest covered by the PDB structure and
        # the ID of the hit of HHsearch corresponding
        # to the current row of the dataframe
        pdb_id = row[PDB_ID_FIELD]
        pdb_chain = row[PDB_CHAIN_FIELD]
        uniprot_segments = row[UNP_SEGMENTS_FIELD]

        if source == "homologs":
            hit_id = index

        # If data for the current PDB are not available
        # yet (it can happen to encounter the same PDB
        # ID because of multiple chains than can
        # represent the protein of interest)
        if not pdb_data[pdb_id]:

            # Try to fetch the PDB data. Use only the first model
            # for multi-model structures
            try:
                pdb_data = pdb_tools.get_pdb_data(pdb_id = pdb_id,
                                                  model = 1)
            
            # If something went wrong, log it and go the next structure
            except Exception as e:
                warnstr = \
                    f"Could not fetch data for the PDB structure " \
                    f"with ID {pdb_id}. This structure will be " \
                    f"ignored. Exception: {e}"
                # Log a warning
                logger.warning(warnstr)
                # Go to the next structure
                continue

            # Add PDB data to the dictionary of PDB data
            pdb_data_dict[pdb_id] = pdb_data
        
        # For each SLiM still not found in a structure
        for slim in still_to_find:

            # Unpack the SLiM data
            s_type, c_start, c_seq, c_end, e_start, e_seq, e_end = \
                slim

            # If it is only the core of the SLiM that needs to
            # be in the PDB structure
            if find_in_pdb == "core":
                # The sequence range will be the one spanning
                # the core motif
                start, seq, end = c_start, c_seq, c_end

            # If it is the extended SLiM that needs to be in
            # the PDB structure
            elif find_in_pdb == "extended":
                # The sequence range will be the one spanning
                # the extended motif
                start, seq, end = e_start, e_seq, e_end
            
            # Try to find out if the SLiM is in the protein
            # segments covered by the PDB structure
            try:
                slim_in_uniprot_segments = \
                    region_in_protein_segments(\
                        seq_range = (start, end),
                        uniprot_segments = uniprot_segments)

            # If something went wrong, log it and go to the next
            # SLiM
            except Exception as e:
                warnstr = \
                    f"Could not check if {start}-{seq}-{end} is in " \
                    f"the segments of the UniProt sequence covered " \
                    f"in the PDB structure with ID {pdb_id}. " \
                    f"This SLiM will not be considered. " \
                    f"Exception: {e}"
                # Log a warning
                logger.warning(warnstr)
                # Go to the next SLiM
                continue

            # If the SLiM is in the structure
            if slim_in_uniprot_segments:

                if source == "protein":
                    
                    # No need to check the PDB file for missing
                    # residues, since all segments with problems
                    # in the mapping have already been excluded
                    
                    # Remove the SLiM from the list of those
                    # to be yet found
                    still_to_find.remove(slim)

                    # Add the SLiM to the list of SLiMs found in
                    # the current PDB structure
                    slims_in_pdb.append([s_type, \
                                         (start, seq, end)])

                elif source == "homologs":
                
                    # Try to find out if the sequence in the homolog
                    # corresponding to the SLiM in the protein
                    # of interest is in the PDB structure
                    try:
                        slim_in_pdb_segments, start_in_t, \
                            seq_in_t, end_in_t = \
                            region_in_homolog_segments(\
                                pdb_data = pdb_data,
                                pdb_id = pdb_id,
                                pdb_chain = pdb_chain,
                                hit_al = als[hit_id],
                                seq_range = (start, end))

                    # If something went wrong, log it and go to the
                    # next SLiM
                    except Exception as e:
                        warnstr = \
                            f"Could not check whether the sequence " \
                            f"{start}-{seq}-{end} is actually in " \
                            f"the PDB structure with ID {pdb_id}. " \
                            f"This SLiM will be ignored. " \
                            f"Exception: {e}"
                        # Log a warning
                        logger.warning(warnstr)
                        # Go to the next SLiM
                        continue

                    # If it is in the PDB structure
                    if slim_in_pdb_segments:
                    
                        # Remove the SLiM from the list
                        # of those to be yet found
                        still_to_find.remove(slim)

                        # Add the SLiM and the sequence corresponding
                        # to the SLiM in the homolog to the list of
                        # SLiMs found in the current PDB structure
                        slims_in_pdb.append(\
                            [s_type, \
                             (start_in_t, seq_in_t, end_in_t), \
                             (start, seq, end)])

        # If no SLiMs are found in the PDB structure
        if not slims_in_pdb:
            # Go on to the next PDB structure
            continue

        # Write out the PDB file corresponding to the structure
        with open(pdb_file, "w") as out:
            out.write(pdb_data)

        # Add useful data to the list containing data about
        # the PDB structure and the SLiMs found in it
        structures_and_slims.append(source, 
                                    pdb_id,
                                    pdb_chain,
                                    pdb_file,
                                    slims_in_pdb) 

    # Return the list of PDB structures and the associated SLiMs
    return structures_and_slims, still_to_find



########################## LOOK FOR DISORDER ##########################



def get_dssp_data(pdb_file,
                  dssp_params,
                  dssp_mapping,
                  wd):
    """Get data from the DSSP output.
    """


    #--------------------------- Run DSSP ----------------------------#


    # Try to run DSSP
    try:
        _, dssp_output = \
            run_tools.run_dssp(pdb_file = pdb_file,
                               dssp_params = dssp_params,
                               wd = wd)

    # If something went wrong, log it and raise an exception
    except Exception as e:
        errstr = \
            f"Could not complete the DSSP run for the PDB file " \
            f"{pdb_file}. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)
    


    #--------------------- Parse the DSSP output ---------------------#


    # Create the DSSP output parser
    parser = dssp_tools.DSSPParser(dssp_mapping = dssp_mapping)

    # Try to parse the DSSP output
    try:
        dssp_data = parser.parse(dssp_output = dssp_output)

    # If something went wrong, log it and raise an exception
    except Exception as e:
        errstr = \
            f"Could not parse the DSSP output file for the PDB " \
            f"file {pdb_file}. Exception: {e}"
        # Log the error
        logger.error(errstr)
        # Re-raise the exception
        raise Exception(errstr)

    # Return the parsed data
    return dssp_data


def is_region_disordered(seq,
                         pdb_id,
                         pdb_chain,
                         fasta_dict,
                         seqdssp_dict,
                         min_len_disorder,
                         disorder_symbols):
    """Check if a portion of the sequence is in a disordered region.
    """
    
    # Get the FASTA sequence of the chain of interest
    chain_f = fasta_dict[chain]
    # Get the FASTA-like DSSP sequence of the chain of interest
    chain_d = seqdssp_dict[chain]

    # Get all the occurrences of the sequence of interest in the
    # FASTA sequence of the chain
    occur = \
        [(i, chain_f[i:i+len(seq)], chain_f[i:i+len(seq)].find(seq)) \
         for i in range(len(chain_f)) \
         if chain_f[i:i+len(seq)].find(seq) != -1]

    # If no occurrence was found, warn the user and return
    if len(occur) == 0:
        warnstr = \
            f"Sequence {seq} not found in PDB {pdb_id} " \
            f"chain {pdb_chain}."
        log.warning(warnstr)
        return
    
    # If multiple occurrences of the sequence were found, return
    # since we have no way of knowing which one corresponds to
    # the one we are interested in. Warn the user.
    elif len(occur) > 1:
        warnstr = \
            f"Multiple occurrences of {seq} found in PDB {pdb_id} " \
            f"chain {pdb_chain} and no way of determining which " \
            f"one is the one of interest. Please check manually."
        log.warning(warnstr)
        return

    # Get the starting position of the only occurrence found
    start_ix = occur[0][2]
    # Get the ending position of the only occurrence found
    end_ix = start_ix+len(seq)

    # Get the FASTA-like DSSP sequence of the occurrence
    seq_seqddsp = chain_d[start_ix:end_ix]

    # We want to be sure that the occurrence is found in
    # a disordered region at least min_len_disorder residues
    # long

    # Subtract the length of the sequence of interest to
    # the minimum length of the disordered region in order
    # to obtain the size of the extension of the sequence
    # on both sides
    side_ext = math.ceil((min_len_disorder-len(seq))/2)
    # If the result is a negative extension (if the minimum
    # length of the disordered region is lower than the
    # length of the sequence of interest)
    if side_ext < 0:
        # The extension will be of length 0
        side_ext = 0

    # Get the putative starting and ending positions of the
    # extended sequence of interest
    ext_start_ix = start_ix - side_ext
    ext_end_ix = end_ix + side_ext

    # Try to extend the sequence on the left side
    try:
        ext_seq_seqddsp = chain_d[ext_start_ix]
    # If there are not enough residues to extend the sequence
    # on the left side by the required number of residues,
    # just extend it as much as possible
    except IndexError:
        warnstr = \
            f"The sequence {seq} is close to the N-terminus of " \
            f"PDB {pdb_id} chain {pdb_chain}, therefore the " \
            f"context on the left side is shorter than " \
            f"{min_len_disorder} residues. The available context " \
            f"will be considered for disorder content."
        log.warning(warnstr)
        # The new starting point will be the N-terminus
        ext_start_ix = 0

    # Try to extend the sequence on the right side
    try:
        ext_seq_seqddsp = chain_d[ext_end_ix]
    # If there are not enough residues to extend the sequence
    # on the right side by the required number of residues,
    # just extend it as much as possible
    except IndexError:
        warnstr = \
            f"The sequence {seq} is close to the C-terminus of " \
            f"PDB {pdb_id} chain {pdb_chain}, therefore the " \
            f"context on the right side is shorter than " \
            f"{min_len_disorder} residues. The available context " \
            f"will be considered for disorder content."
        log.warning(warnstr)
        # The new ending point will be the N-terminus
        ext_end_ix = len(chain_d) - 1

    # Get the extended sequence
    ext_seq = chain_d[ext_start_ix:ext_end_ix]

    # If all residues in the extended sequence are disordered
    if set(ext_seq).issubset(set(disorder_symbols)):
        return True
    # Otherwise
    else:
        return False


def get_results(structures_and_slims,
                wd):
    """Get the results for the SLiMs found."""

    # Create a dictionary to store the SLiMs found and annotated
    slims_results = {"protein" : [], "homologs" : []}

    for source, pdb_id, pdb_chain, pdb_file, slims_in_pdb \
        in structures_and_slims:

        # Try to run DSSP, parse the DSSP output and get the data
        try:
            _, fasta_dict, seqdssp_dict = \
                get_dssp_data(pdb_file = pdb_file,
                              dssp_params = dssp_params,
                              dssp_mapping = dssp_mapping,
                              wd = wd)

        # If something went wrong, log it and go to the next structure
        except Exception as e:
            warnstr = \
                f"Could not get DSSP data for the PDB structure " \
                f"with ID {pdb_id}. This structure will be " \
                f"ignored. Exception: {e}"
            # Log a warning
            logger.warning(warnstr)
            # Go to the next structure
            continue

        # For each SLiM found in the PDB structure
        for slim_data in slims_in_pdb:
    
            # If the SLiM comes from a structure of the protein
            # of interest
            if source == "protein":
                # Data will only be about the SLiM in the protein
                slim_type, (start, seq, end) = slim_data
                start_q, seq_q, end_q = start, seq, end
                
            # If the SLiM comes from a structure of a homolog
            elif source == "homologs":
                # Data will be both about the SLiM in the protein
                # and about its corresponding sequence in the
                # homolog ("_q" stands for "query", therefore the
                # protein of interest)
                slim_type, (start, seq, end), \
                    (start_q, seq_q, end_q) = slim_data

            # Try to check if the SLiM is in a disordered region
            try:
                is_disordered = \
                    is_region_disordered(\
                        seq = seq,
                        pdb_id = pdb_id,
                        pdb_chain = pdb_chain,
                        fasta_dict = fasta_dict,
                        seqdssp_dict = seqdssp_dict,
                        min_len_disorder = min_len_disorder,
                        disorder_symbols = disorder_symbols)

            # If something went wrong, log it and go to the 
            # next structure
            except Exception as e:
                warnstr = \
                    f"Could not check whether the SLiM with Uniprot " \
                    f"sequence {start_q}-{seq_q}-{end_q} is in a " \
                    f"disordered region in the PDB structure with " \
                    f"ID {pdb_id} ({source}), chain {pdb_chain}. " \
                    f"This SLiM will be ignored. Exception: {e}"
                # Log a warning
                logger.warning(warnstr)
                # Go to the next structure
                continue

            # If the SLiM comes from a structure of the protein
            # of interest
            if source == "protein":
                # Save only information about the SLiM in the protein
                results = {"slim_type" : slim_type,
                           "slim_start" : start_q,
                           "slim_sequence" : seq_q, 
                           "slim_end" : end_q,
                           "disordered_region" : str(is_disordered)}

            # If the SLiM comes from a structure of a homolog
            elif source == "homologs":
                # Save information about both the SLiM in the
                # protein and its corresponding sequence in the
                # homolog
                results = {"slim_type" : slim_type,
                           "slim_start" : start_q,
                           "slim_sequence" : seq_q, 
                           "slim_end" : end_q,
                           "slim_start_pdb" : start,
                           "slim_sequence_pdb" : seq,
                           "slim_end_pdb" : end,
                           "disordered_region" : str(is_disordered)}

            # Add the SLiM to results
            slims_results[source].append(results)

    # Get the results' dataframes
    results_protein_df = pd.DataFrame(slims_results["protein"])
    results_homologs_df = pd.DataFrame(slims_results["homologs"])

    # Return the results' dataframes
    return results_protein_df, results_homologs_df



def get_disprot_annotations(uniprot_id,
                            slims_to_find):
    """Get the DisProt annotations for a set of SLiMs, given
    the UniProt ID to query for them.
    """

    # Create the parser for DisProt annotations
    parser = disprot_tools.DisProtParser()

    # Get the annotations
    df = parser.get_annotations(uniprot_id = uniprot_id)

    # Get only the portion of the dataframe regarding disordered
    # regions, and only select the columns containing the starting
    # and ending points of the regions
    disorder_df = \
        df.loc[df["term_name"] == "Disorder"][["start", "end"]]

    # Create an empty dictionary to store the results
    results = {}

    # For each SLiM
    for slim in slims_to_find:

        # Unpack the SLiM data
        s_type, c_start, c_seq, c_end, e_start, e_seq, e_end = slim

        # If it is only the core of the SLiM that needs to
        # be in the PDB structure
        if find_in_pdb == "core":
            # The sequence range will be the one spanning
            # the core motif
            seq_range = (c_start, c_end)
            seq = c_seq

        # If it is the extended SLiM that needs to be in
        # the PDB structure
        elif find_in_pdb == "extended":
            # The sequence range will be the one spanning
            # the extended motif
            seq_range = (e_start, e_end)
            seq = e_seq
            
        # Try to find out if the SLiM is in one of the
        # disordered regions

        # Initialize the disordered flag to False
        is_disordered = False

        # For each disordered region
        for row in disorder_df.iterrows():
            # If the SLiM is included in the disordered region
            if start >= row["start"] and end <= row["end"]:
                is_disordered = True

        # Save information about the SLiM
        results.append({"slim_type" : s_type,
                        "slim_start" : seq_range[0],
                        "slim_sequence" : seq, 
                        "slim_end" : seq_range[1],
                        "disordered_region" : str(is_disordered)})

    # Convert the results into a daraframe and return them
    return pd.DataFrame(results)

    
