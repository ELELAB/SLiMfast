#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    util.py
#
#    Miscellanea utility functions for the SLiMfast pipeline.
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
import math
import os
import re
import requests
import sys
# Third-party packages
import pandas as pd
import yaml
# SLiMfast
from .defaults import (
    ALPHAFOLD_ID_FIELD,
    CONFIG_RUN_DIR,
    DEPOSITION_DATE_FIELD,
    HIT_ID_FIELD,
    METHOD_FIELD,
    RANGE_SEP,
    RESOLUTION_FIELD,
    SEGMENTS_SEP,
    STRUCT_CHAIN_FIELD,
    STRUCT_FRAGMENT_FIELD,
    STRUCT_ID_FIELD,
    STRUCT_SEGMENTS_FIELD,
    STRUCT_VERSION_FIELD,
    UNP_ID_FIELD,
    UNP_ISOFORM_FIELD,
    UNP_SEGMENTS_FIELD)
from . import (
    alphafold_tools,
    dataframe_tools,
    disprot_tools,
    dssp_tools,
    hhsearch_tools,
    pdb_tools,
    uniprot_tools,
    run_tools)



# Get the module logger
logger = log.getLogger(__name__)



########################## GENERAL UTILITIES ##########################



def get_wd(config_run):
    """Get the directory where the pipeline will be run.
    """

    # Get the working directory from the configuration file
    wd = config_run["wd"]
    
    # If the working directory is defined as an empty string
    if wd == "":
        
        # The working directory will be the current working
        # directory
        wd = os.getcwd()
    
    else:
        
        # If it is a directory name
        if os.path.basename(wd) == wd:
            
            # The working directory will be a new directory
            # with that name inside the current working
            # directory
            wd = os.path.join(os.getcwd(), wd)
        
        # If it is a directory path
        else:
            
            # Get the absolute path of the directory and set
            # it as the working directory
            wd = os.path.abspath(wd)

    # Return the working directory
    return wd


def get_abspath(path):
    """Given a path, return its absolute path. Return
    None if the path given is None.
    """

    return os.path.abspath(path) if path is not None else path


def get_items(d,
              keys,
              expected_num_values,
              insert_first = None):
    """Get multiple items from a dictionary, given their keys,
    and return a list of them.
    """

    # Create an empty list of itens
    items = []

    # If there are elements to be inserted first
    if insert_first is not None:

        # Add them
        items.extend(insert_first)

    # For each key
    for k, num_v in zip(keys, expected_num_values):

        # If the key is in the dictionary
        if k in d.keys():
            
            # Get the corresponding value
            v = d[k]

            # If the value is a list or a tuple
            if isinstance(v, (list, tuple)):

                # Extend the list
                items.extend(v)

            # Otherwise
            else:

                # Append the value to the list
                items.append(v)

        # Otherwise
        else:

            # Append None to the list for as many
            # times as the number of expected values
            items.extend([None]*num_v)

    # Return the list
    return items



############################## READ FILES #############################



def read_fasta(fasta_file):
    """Read the header and sequence from a FASTA file.
    """

    # Open the FASTA file
    with open(fasta_file, "r") as f:

        # Initialize the header to an empty string
        header = ""
        
        # Initialize the sequence to an empty string
        sequence = ""
        
        # For each line
        for line in f:
            
            # If you encounter an empty line or a comment line
            if re.match(r"^\s*$", line) or line.startswith("#"):
                
                # Ignore it
                continue

            # If the line is the header of the FASTA file
            elif line.startswith(">"):
                
                # Get the header (with the initial '>' stripped)
                header = line.lstrip(">").rstrip("\n")

            # If the line is any other line (sequence line)
            else:
                
                # Add the sequence on the line to the full sequence
                # (the sequence may be split on multiple lines)
                sequence += line.rstrip("\n")

        # Return header and sequence
        return header, sequence


def read_slims_file(slims_file):
    """Get the instances of a SLiM in a sequence from a text file.
    """

    # Open the file containing the SLiM instances
    with open(slims_file, "r") as f:

        # Create an empty list to store the SLiM instances
        slim_instances = []

        # For each line
        for line in f:
            
            # If you encounter an empty line or a comment line
            if re.match(r"^\s*$", line) or line.startswith("#"):

                # Ignore it
                continue
            
            # Unpack the components of each SLiM instance
            # and add the instance to the list
            else:
                s_type, s_start, s_seq, s_end = \
                    line.rstrip("\n").split(",")
                slim_instances.append((s_type, int(s_start),
                                       s_seq, int(s_end)))

        # Return the list of SLiM instances
        return slim_instances


def read_config_file(config_file):
    """Get the configuration from the configuration file.
    """

    # Get the name of the configuration file without extension
    # and the extension
    config_name, config_ext = \
        os.path.splitext(os.path.basename(config_file))

    # If the configuration file is a name without extension
    if config_file == config_name:

        # Assume it is a file with that name inside the default
        # directory for configuration files
        config_file = \
            os.path.join(CONFIG_RUN_DIR, config_name + config_ext)

    # Otherwise
    else:

        # Return the absolute path to the configuration file
        config_file = os.path.abspath(config_file)

    # Load the configuration file
    return yaml.safe_load(open(config_file, "r"))


def read_ids_slims_file(ids_slims_file):
    """Get the UniProt IDs of the proteins of interest and
    associated SLiMs.
    """

    # Read in the CSV file as a data frame
    df = pd.read_csv(ids_slims_file,
                     sep = ",")

    # If there are less than the expected five columns
    # in the data frame
    if len(df.columns) < 5:

        # If there is only one column and it contains the
        # UniProt IDs (which can happen when the SLiM search
        # is not requested as part of the pipeline)
        if len(df.columns) == 1 \
        and UNP_ID_FIELD in df.columns.tolist():

            # Return a data frame mapping the UniProt IDs
            # to None values
            return {k : None for k in df[UNP_ID_FIELD]}

        # Otherwise
        else:

            # Inform the user that there is something wrong
            # with the file and raise an exception
            errstr = \
                f"The file {ids_slims_file} contains less than " \
                f"the expected five columns. The only case in which " \
                f"this is accepted is when the data frame only " \
                f"contains the column '{UNP_ID_FIELD}' with the " \
                f"UniProt IDs of the proteins of interest."
            raise ValueError(errstr)

    # Convert the data frame into a dictionary keyed on the
    # unique UniProt IDs (SLiMs for the same IDs are, therefore,
    # grouped)
    d = df.groupby(UNP_ID_FIELD).agg(\
            lambda x: x.tolist()).to_dict("index")

    # Return a dictionary where keys are the UniProt IDs and values
    # are data frames containing the SLiMs for each UniProt ID
    return {k : pd.DataFrame(v) for k, v in d.items()}



############################ GET STRUCTURES ###########################



def get_structures_df(uniprot_id):
    """Create a data frame containing information about the
    structures available for a protein in the Protein Data Bank,
    given the protein's UniProt ID.
    """
    
    # Create an empty list where to store the data that
    # will be converted into a data frame
    data = []

    # Get the PDB IDs corresponding to the given UniProt ID
    pdbs = uniprot_tools.get_pdbs(uniprot_id = uniprot_id)
    
    # For each PDB ID corresponding to the given UniProt ID
    for struct_id in pdbs:
        
        # Get the mapping between the segments of the protein
        # found in the PDB structure and their corresponding
        # UniProt numbering
        mappings = \
            pdb_tools.get_protein_segments(struct_id = struct_id,
                                           uniprot_id = uniprot_id)
        
        # For each chain and the corresponding UniProt segments
        for struct_chain, segments in mappings.items():
            
            # Get the starting and ending points of the segments of
            # the proteins that are covered, both in the UniProt
            # sequence and in the structure's sequence
            unp_segments = \
                [RANGE_SEP.join([str(t) for t in s[1]]) \
                 for s in segments]
            struct_segments = \
                [RANGE_SEP.join([str(t) for t in s[0]]) \
                 for s in segments]

            # Join the segments in a formatted string
            unp_segments_joined = SEGMENTS_SEP.join(unp_segments)
            struct_segments_joined = SEGMENTS_SEP.join(struct_segments)
            
            # Append the data of interest to the dictionary
            data.append({STRUCT_ID_FIELD : struct_id,
                         STRUCT_CHAIN_FIELD : struct_chain,
                         STRUCT_SEGMENTS_FIELD : struct_segments_joined,
                         UNP_ID_FIELD : uniprot_id,
                         UNP_SEGMENTS_FIELD : unp_segments_joined})

    # If no matching PDB structures were found for the
    # UniProt ID
    if not data:
        
        # Return an empty data frame
        return pd.DataFrame()

    # Create a data frame from the data
    df = pd.DataFrame(data)

    # Return the data frame
    return df


def get_alphafold_df(uniprot_id):
    """Create a data frame containing information about the
    structures available for a protein in the AlphaFold Protein
    Structure Database, given the protein's UniProt ID.
    """
    
    # Create an empty list where to store the data that
    # will be converted into a data frame
    data = []

    # Get the PDB IDs corresponding to the given UniProt ID
    structs = \
        alphafold_tools.get_alphafold_data(uniprot_id = uniprot_id)

    # If no matching PDB structures were found for the
    # UniProt ID
    if not structs:
        
        # Return an empty data frame
        return pd.DataFrame()
    
    # For each structure ID corresponding to the given UniProt ID
    for struct_id, (version, fragment, _) in structs.items():
            
        # Append the data of interest to the dictionary
        data.append({STRUCT_ID_FIELD : struct_id,
                     UNP_ID_FIELD : uniprot_id,
                     STRUCT_FRAGMENT_FIELD : fragment,
                     STRUCT_VERSION_FIELD : version})

    # Create a data frame from the data
    df = pd.DataFrame(data)

    # Return the data frame
    return df


def add_structures_metadata(df,
                            source):
    """Add metadata regarding each structure to the
    data frame.
    """

    # Create an empty list to store the metadata
    # for each unique structure ID
    metadata_list = []

    # For each unique ID
    for struct_id in df[STRUCT_ID_FIELD].unique():

        # If the structure comes from the Protein Data Bank
        if source == "pdb":

            # Try to get the structure metadata
            try:
                
                metadata = \
                    pdb_tools.get_structure_metadata(\
                        struct_id = struct_id)

            # If something went wrong       
            except Exception as e:

                # Log a warning
                warnstr = \
                    f"Could not fetch metadata for the structure " \
                    f"with ID {struct_id} from the Protein Data " \
                    f"Bank. This structure will be ignored. " \
                    f"Exception: {e}"
                logger.warning(warnstr)
                
                # Go to the next structure
                continue

        # If the structure comes from the AlphaFoldDB
        elif source == "alphafold":
        
            # Try to get the structure metadata
            try:
                
                metadata = \
                    alphafold_tools.get_structure_metadata(\
                        struct_id = struct_id)

            # If something went wrong       
            except Exception as e:

                # Log a warning
                warnstr = \
                    f"Could not fetch metadata for the structure " \
                    f"with name {struct_id} from the AlphaFold " \
                    f"Protein Structure Database. This structure " \
                    f"will be ignored. Exception: {e}"
                logger.warning(warnstr)
                
                # Go to the next structure
                continue


        # Add a field for the structure ID (it is needed to merge
        # the dataframes later)
        metadata.update({STRUCT_ID_FIELD : struct_id})

        # Add the dictionary of metadata to the list
        metadata_list.append(metadata)

        # Inform the user that the metadata were collected
        # successfully for the current structure
        infostr = \
            f"Metadata successfully collected for the structure " \
            f"with ID {struct_id}."
        logger.info(infostr)


    # If metadata for at least one of the structures were retrieved
    if metadata_list:

        # Create a data frame containing the metadata for all
        # structures
        metadata_df = pd.DataFrame(metadata_list)

        # Combine the concatenated data frame with the data frame
        # containing the structures' metadata
        updated_df = pd.merge(left = df,
                              right = metadata_df,
                              on = STRUCT_ID_FIELD,
                              how = "outer")

        # Return the updated data frame
        return updated_df


    # Otherwise
    else:

        # Inform the user that no metadata were found for the
        # structures
        warnstr = \
            "No metadata could be added for any structure."
        logger.warning(warnstr)

        # Return the original data frame
        return df



############################## FIND SLIMS #############################



def region_in_protein_segments(seq_range,
                               uniprot_segments,
                               struct_segments):
    """Check if a sequence belongs at least to one UniProt
    segment, given the sequence numeric range.
    """
    
    # Get the starting and ending point of the sequence range  
    start_s, end_s = seq_range

    # Get the length of the sequence
    length = len(range(start_s, end_s))

    # Try to get the numeric ranges of all UniProt segments passed
    try:
        
        uniprot_ranges = \
            [(int(s.split(RANGE_SEP)[0]), int(s.split(RANGE_SEP)[1])) \
             for s in uniprot_segments.split(SEGMENTS_SEP)]
    
    # If something went wrong
    except ValueError:

        # Log the error
        errstr = \
            "Could not find the UniProt segments covered by the " \
            "structure, since at least one extreme of at least " \
            "one of the intervals defining the segments is None."

        # Re-raise the exception
        raise ValueError(errstr)

    # Try to get the numeric ranges of all PDB segments passed
    try:
        
        struct_ranges = \
            [(int(s.split(RANGE_SEP)[0]), int(s.split(RANGE_SEP)[1])) \
             for s in struct_segments.split(SEGMENTS_SEP)]
    
    # If something went wrong
    except ValueError:

        # Log the error
        errstr = \
            "Could not find the segments covered by the " \
            "structure, since at least one extreme of at least " \
            "one of the intervals defining the segments is None."
        
        # Re-raise the exception
        raise ValueError(errstr)

    # For each UniProt segment
    for uniprot_start, uniprot_end in uniprot_ranges:

        # If the SLiM is found within the segment
        if start_s >= uniprot_start and end_s <= uniprot_end:
            
            # Return True                   
            return True

    # If the sequence was never found, return False
    return False


def region_in_homolog_segments(pdb_data,
                               struct_id,
                               struct_chain,
                               hit_al,
                               seq_range):
    """Check if a region of the protein of interest is found in the
    segments corresponding to the sequences aligned to such a
    region in the proteins' homolog.
    """
    

    #-------------------- Find the aligned sequence ------------------#


    # Try to find the sequence in the homolog aligned to the
    # sequence of interest in the HHsearch alignments
    try:
        
        start_in_t, slim_in_t, end_in_t, is_match = \
            hhsearch_tools.get_aligned_sequence(hit_al = hit_al,
                                                seq_range = seq_range)
    
    # If something went wrong
    except Exception as e:

        # Log the error
        start, end = seq_range
        errstr = \
            f"Could not find the sequence {start}-{end} in the " \
            f"alignment. Exception: {e}"
        
        # Re-raise the exception
        raise Exception(errstr)

    # If the aligned sequence does not perfectly match the sequence
    # of interest
    if not is_match:

        # Return
        return False, None, None, None

    # Get the numeric range covered by the aligned sequence in
    # the PDB structure of the homolog
    slim_range = set(list(range(start_in_t, end_in_t)))


    #------------------- Find the missing residues -------------------#


    # Try to get the missing residues in the structure
    try:
        missing_res = \
            pdb_tools.get_missing_residues(pdb_data = pdb_data)
    
    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not get the missing residues for the PDB " \
            f"structure with ID {struct_id}. Exception: {e}"
        
        # Re-raise the exception
        raise Exception(errstr)
    
    # Get the missing residues for the chain of interest
    missing_res_chain = set(missing_res[struct_chain])

    # If none of the residues in the aligned sequence of interest
    # is missing in the structure of the homolog
    if not set(slim_range).intersection(missing_res_chain):

        # Return that none of the residues were missing,
        # plus the start, sequence, and end of the SLiM
        # in the region
        return True, start_in_t, slim_in_t, end_in_t

    # Return None for each of the structure-related info if
    # some residues of the aligned sequence of interest were
    # missing in the structure of the homolog
    return True, None, None, None


def find_slims(source,
               structures_df,
               slims_to_find,
               struct_dir,
               als = None):
    """Find the SLiMs in the structures of the protein
    of interest or in those of the homologs.
    """
        
    # Create an empty dictionary storing the data
    # for the structures of interest
    struct_data_dict = {}

    # Create an empty list to store the PDB files
    # and the SLiMs to be checked for disorder
    structures_and_slims = []

    # Create a set to store the SLiMs found
    slims_found = set()

    # For each row in the data frame (= structure)
    for index, row in structures_df.iterrows():

        # If no SLiMs are left to search
        if set(slims_to_find) == set(slims_found):

            # Exit the loop
            break

        # Create an empty list to store the SLiMs
        # found in each structure
        slims_in_struct = []

        # Get the structure ID, chain, the segments of the protein
        # of interest covered by the structure
        struct_id = row[STRUCT_ID_FIELD]
        struct_chain = row[STRUCT_CHAIN_FIELD]
        uniprot_segments = row[UNP_SEGMENTS_FIELD]
        struct_segments = row[STRUCT_SEGMENTS_FIELD]

        # If the structure is a homolog
        if source == "homologs":

            # The hit IDs will be the index of the data frame
            hit_id = index

        # If data for the current structure are not available yet (it
        # can happen to encounter the same ID because of multiple
        # chains than can represent the protein of interest)
        if not struct_id in struct_data_dict.keys():

            # If the structure comes from the Protein Data Bank
            if source in ("structures", "homologs"):

                # Try to fetch the structure's data. Use only the first
                # model for multi-model structures
                try:
                    
                    struct_data = \
                        pdb_tools.get_pdb_data(struct_id = struct_id,
                                               model = 1)

                    # Inform the user that data for the current structure
                    # were successfully retrieved
                    infostr = \
                        f"Data for the structure with ID {struct_id} " \
                        f"have been successfully retrieved for the " \
                        f"Protein Data Bank."
                    logger.info(infostr)
                
                # If something went wrong
                except Exception as e:

                    # Log a warning
                    warnstr = \
                        f"Could not fetch data for the structure " \
                        f"with ID {struct_id} from the Protein Data " \
                        f"Bank. This structure will be ignored. " \
                        f"Exception: {e}"
                    logger.warning(warnstr)
                    
                    # Go to the next structure
                    continue

            # If the structure comes from the AlphaFold Protein
            # Structure Database
            elif source == "alphafold":

                # Try to fetch the structure's data.
                try:
                    
                    struct_data = \
                        alphafold_tools.get_alphafold_data(\
                            struct_id = struct_id)

                    # Inform the user that data for the current structure
                    # were successfully retrieved
                    infostr = \
                        f"Data for the structure with ID {struct_id} " \
                        f"have been successfully retrieved for the " \
                        f"AlphaFold Protein Structure Database."
                    logger.info(infostr)
                
                # If something went wrong
                except Exception as e:

                    # Log a warning
                    warnstr = \
                        f"Could not fetch data for the structure " \
                        f"with ID {struct_id} from the AlphaFold " \
                        f"Protein Structure Database. This " \
                        f"structure will be ignored. " \
                        f"Exception: {e}"
                    logger.warning(warnstr)
                    
                    # Go to the next structure
                    continue

            # Add data to the dictionary of structures' data
            struct_data_dict[struct_id] = struct_data
        
        # For each SLiM still not found in a structure
        for slim in slims_to_find:

            # Ignore the SLiM if it has been found already
            if slim in slims_found:

                # Go to the next SLiM
                continue

            # Unpack the SLiM data
            s_type, s_start, s_seq, s_end = slim
            
            # Try to find out if the SLiM is in the protein
            # segments covered by the structure
            try:
                
                slim_in_uniprot_segments = \
                    region_in_protein_segments(\
                        seq_range = (s_start, s_end),
                        uniprot_segments = uniprot_segments,
                        struct_segments = struct_segments)

            # If something went wrong
            except Exception as e:

                # Log a warning
                warnstr = \
                    f"Could not check if {s_start}-{s_seq}-{s_end} " \
                    f"is in the segments of the UniProt sequence " \
                    f"covered in the structure with ID {struct_id}. " \
                    f"This SLiM will be ignored for this structure."
                logger.warning(warnstr)
                
                # Go to the next SLiM
                continue

            # If the SLiM is in the structure
            if slim_in_uniprot_segments:

                # If the structure is of the protein of interest
                if source in ("structures", "alphafold"):
                    
                    # No need to check the structure for missing
                    # residues, since all segments with problems
                    # in the mapping have already been excluded
                    
                    # Add the SLiM to the set of SLiMs found
                    slims_found.add(slim)

                    # Add the SLiM to the list of SLiMs found in
                    # the current structure
                    slims_in_struct.append([s_type, \
                                           (s_start, s_seq, s_end)])

                # If the structure is of a homolog
                elif source == "homologs":
                
                    # Try to find out if the sequence in the homolog
                    # corresponding to the SLiM in the protein
                    # of interest is in the structure
                    try:
                        slim_in_struct_segments, start_in_t, \
                            seq_in_t, end_in_t = \
                            region_in_homolog_segments(\
                                pdb_data = pdb_data,
                                struct_id = struct_id,
                                struct_chain = struct_chain,
                                hit_al = als[hit_id],
                                seq_range = (start, end))

                    # If something went wrong
                    except Exception as e:

                        # Log a warning
                        warnstr = \
                            f"Could not check whether the sequence " \
                            f"{s_start}-{s_seq}-{s_end} is actually " \
                            f"in the PDB structure with ID " \
                            f"{struct_id}. This SLiM will be ignored " \
                            f"for this structure."
                        logger.warning(warnstr)
                        
                        # Go to the next SLiM
                        continue

                    # If the SLiM is in the structure
                    if slim_in_struct_segments:
                    
                        # Add the SLiM to the set of SLiMs found
                        slims_found.add(slim)

                        # Add the SLiM and the sequence corresponding
                        # to the SLiM in the homolog to the list of
                        # SLiMs found in the current structure
                        slims_in_struct.append(\
                            [s_type, \
                             (start_in_t, seq_in_t, end_in_t), \
                             (s_start, s_seq, s_end)])

        # If no SLiMs are found in the structure
        if not slims_in_struct:

            # Inform the user of it
            infostr = \
                f"No SLiMs found in the structure with ID {struct_id}."
            logger.info(infostr)
            
            # Go on to the next structure
            continue

        # Inform the user of the SLiMs found in the current
        # structure
        slims_str = \
            ", ".join(\
                [f"{s_type}_{s_start}-{s_seq}-{s_end}" \
                 for (s_type, (s_start, s_seq, s_end)) \
                 in slims_in_struct])
        infostr = \
            f"The following SLiMs have been found in the structure " \
            f"with ID {struct_id}: {slims_str}."
        logger.info(infostr)

        # Write out the file corresponding to the structure (it
        # is written out only if at least a SLiM has been found
        # in it)
        os.makedirs(struct_dir, exist_ok = True)
        struct_file = os.path.join(struct_dir, f"{struct_id}.pdb")
        with open(struct_file, "w") as out:
            out.write("\n".join(struct_data))

        # Inform the user that the structure was saved successfully
        infostr = \
            f"The structure with ID {struct_id} has been " \
            f"written successfully to {struct_file}."
        logger.info(infostr)

        # Add useful data to the list containing data about
        # the structure and the SLiMs found in it
        structures_and_slims.append((source, 
                                     struct_id,
                                     struct_chain,
                                     struct_file,
                                     slims_in_struct))

    # Return the list of PDB structures and the associated SLiMs
    return structures_and_slims, \
           [s for s in slims_to_find if s not in slims_found]



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

        # Inform the user that the DSSP run was successful
        infostr = \
            f"DSSP was successfully run on the PDB file {pdb_file}."

    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not complete the DSSP run for the PDB file " \
            f"{pdb_file}. Exception: {e}"
        logger.exception(errstr)
        
        # Re-raise the exception
        raise Exception(errstr)
    


    #--------------------- Parse the DSSP output ---------------------#


    # Create the DSSP output parser
    parser = dssp_tools.DSSPParser(dssp_mapping = dssp_mapping)

    # Try to parse the DSSP output
    try:
        
        dssp_data = parser.parse(dssp_output = dssp_output)

        # Inform the user of the successful creation of the
        # data frame
        infostr = \
            f"Data in the DSSP output {dssp_output} have been " \
            f"parsed successfully."
        logger.info(infostr)

    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not parse the DSSP output file for the PDB " \
            f"file {pdb_file}. Exception: {e}"
        logger.exception(errstr)
        
        # Re-raise the exception
        raise Exception(errstr)

    # Return the parsed data
    return dssp_data


def is_region_disordered(seq,
                         struct_id,
                         struct_chain,
                         fasta_dict,
                         seqdssp_dict,
                         min_len_disorder,
                         disorder_symbols):
    """Check if a portion of the sequence is in a disordered region.
    """
    
    # Get the FASTA sequence of the chain of interest
    chain_f = fasta_dict[struct_chain]
    
    # Get the FASTA-like DSSP sequence of the chain of interest
    chain_d = seqdssp_dict[struct_chain]

    # Get all the occurrences of the sequence of interest in the
    # FASTA sequence of the chain
    occur = \
        [(i, chain_f[i:i+len(seq)], chain_f[i:i+len(seq)].find(seq)) \
         for i in range(len(chain_f)) \
         if chain_f[i:i+len(seq)].find(seq) != -1]

    # If no occurrence was found, warn the user and return
    if len(occur) == 0:
        warnstr = \
            f"Sequence {seq} not found in PDB {struct_id} " \
            f"chain {struct_chain}."
        logger.warning(warnstr)
        return
    
    # If multiple occurrences of the sequence were found, return
    # since we have no way of knowing which one corresponds to
    # the one we are interested in. Warn the user.
    elif len(occur) > 1:
        warnstr = \
            f"Multiple occurrences of {seq} found in structure " \
            f"{struct_id}, chain {struct_chain}, and SLiMfast has " \
            f"no way of determining which one is the one of " \
            f"interest. Please check manually."
        logger.warning(warnstr)
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

        # Inform the user of the extension
        warnstr = \
            f"The sequence {seq} is close to the N-terminus of " \
            f"PDB {struct_id} chain {struct_chain}, therefore the " \
            f"context on the left side is shorter than " \
            f"{min_len_disorder} residues. The available context " \
            f"will be considered for disorder content."
        logger.warning(warnstr)
        
        # The new starting point will be the N-terminus
        ext_start_ix = 0

    # Try to extend the sequence on the right side
    try:
        ext_seq_seqddsp = chain_d[ext_end_ix]
    
    # If there are not enough residues to extend the sequence
    # on the right side by the required number of residues,
    # just extend it as much as possible
    except IndexError:

        # Inform the user of the extension
        warnstr = \
            f"The sequence {seq} is close to the C-terminus of " \
            f"PDB {struct_id} chain {struct_chain}, therefore the " \
            f"context on the right side is shorter than " \
            f"{min_len_disorder} residues. The available context " \
            f"will be considered for disorder content."
        logger.warning(warnstr)
        
        # The new ending point will be the N-terminus
        ext_end_ix = len(chain_d) - 1

    # Get the extended sequence
    ext_seq = chain_d[ext_start_ix:ext_end_ix]

    # If all residues in the extended sequence are disordered
    if set(ext_seq).issubset(set(disorder_symbols)):

        # Return True (the region is disordered)
        return True
    
    # Otherwise
    else:

        # Return False (the region is not disordered)
        return False


def get_results(structures_and_slims,
                dssp_params,
                dssp_mapping,
                min_len_disorder,
                disorder_symbols,
                wd):
    """Get the results for the SLiMs found."""

    # Create a dictionary to store the SLiMs found and annotated
    slims_results = \
        {"structures" : [], "alphafold" : [], "homologs" : []}

    # For each structure and SLiMs found in it
    for source, struct_id, struct_chain, pdb_file, slims_in_struct \
        in structures_and_slims:

        # Try to run DSSP, parse the DSSP output and get the data
        try:
            _, fasta_dict, seqdssp_dict = \
                get_dssp_data(pdb_file = pdb_file,
                              dssp_params = dssp_params,
                              dssp_mapping = dssp_mapping,
                              wd = wd)

        # If something went wrong
        except Exception as e:

            # Log a warning
            warnstr = \
                f"Could not get DSSP data for the structure " \
                f"with ID {struct_id}. This structure will be " \
                f"ignored. Exception: {e}"
            logger.warning(warnstr)
            
            # Go to the next structure
            continue

        # For each SLiM found in the structure
        for slim_data in slims_in_struct:
    
            # If the SLiM comes from a structure of the protein
            # of interest
            if source in ("structures", "alphafold"):
                
                # Data will only be about the SLiM in the protein
                slim_name, (start, seq, end) = slim_data
                start_q, seq_q, end_q = start, seq, end
                
            # If the SLiM comes from a structure of a homolog
            elif source == "homologs":
                
                # Data will be both about the SLiM in the protein
                # and about its corresponding sequence in the
                # homolog ("_q" stands for "query", therefore the
                # protein of interest)
                slim_name, (start, seq, end), \
                    (start_q, seq_q, end_q) = slim_data

            # Try to check if the SLiM is in a disordered region
            try:
                
                is_disordered = \
                    is_region_disordered(\
                        seq = seq,
                        struct_id = struct_id,
                        struct_chain = struct_chain,
                        fasta_dict = fasta_dict,
                        seqdssp_dict = seqdssp_dict,
                        min_len_disorder = min_len_disorder,
                        disorder_symbols = disorder_symbols)

            # If something went wrong
            except Exception as e:

                # Log a warning
                warnstr = \
                    f"Could not check whether the SLiM with UniProt " \
                    f"sequence {start_q}-{seq_q}-{end_q} is in a " \
                    f"disordered region in the structure with " \
                    f"ID {struct_id} ({source}), chain {struct_chain}. " \
                    f"This SLiM will be ignored. Exception: {e}"
                logger.warning(warnstr)
                
                # Go to the next structure
                continue

            # If the SLiM comes from a structure of the protein
            # of interest
            if source in ("structures", "alphafold"):
                
                # Save only information about the SLiM in the protein
                results = {"source" : source,
                           "struct_id" : struct_id,
                           "struct_chain" : struct_chain,
                           "slim_name" : slim_name,
                           "slim_start" : start_q,
                           "slim_sequence" : seq_q, 
                           "slim_end" : end_q,
                           "disordered_region" : str(is_disordered)}

            # If the SLiM comes from a structure of a homolog
            elif source == "homologs":
                
                # Save information about both the SLiM in the
                # protein and its corresponding sequence in the
                # homolog
                results = {"source" : source,
                           "struct_id" : struct_id,
                           "struct_chain" : struct_chain,
                           "slim_name" : slim_name,
                           "slim_start" : start_q,
                           "slim_sequence" : seq_q, 
                           "slim_end" : end_q,
                           "slim_start_pdb" : start,
                           "slim_sequence_pdb" : seq,
                           "slim_end_pdb" : end,
                           "disordered_region" : str(is_disordered)}

            # Add the SLiM to results
            slims_results[source].append(results)

    # Get the results' data frames
    results_protein_df = pd.DataFrame(slims_results["structures"])
    results_alphafold_df = pd.DataFrame(slims_results["alphafold"])
    results_homologs_df = pd.DataFrame(slims_results["homologs"])

    # Return the results' data frames
    return results_protein_df, \
           results_alphafold_df, \
           results_homologs_df


def get_disprot_annotations(uniprot_id,
                            slims_to_find):
    """Get the DisProt annotations for a set of SLiMs, given
    the UniProt ID to query for them.
    """

    # Create the parser for DisProt annotations
    parser = disprot_tools.DisProtParser()

    # Try to get the annotations
    try:

        df = parser.get_annotations(uniprot_id = uniprot_id)

        # Inform the user that the annotations were successfully
        # retrieved
        infostr = \
            f"DisProt annotations for the protein with UniProt " \
            f"ID {uniprot_id} have been successfully retrieved."
        logger.info(infostr)

    # If something went wrong
    except Exception as e:

        # Log it
        errstr = \
            f"Could not get DisProt annotations for the protein " \
            f"with UniProt ID {uniprot_id}. Exception: {e}"
        logger.exception(errstr)

        # Re-raise the exception
        raise Exception(errstr)

    # If there were no annotations
    if df.empty:

        # Return the empty data frame
        return df

    # Get only the portion of the dataframe regarding disordered
    # regions, and only select the columns containing the starting
    # and ending points of the regions
    disorder_df = \
        df.loc[df["term_name"] == "Disorder"][["start", "end"]]

    # Create an empty list to store the results
    results = []

    # For each SLiM
    for slim in slims_to_find:

        # Unpack the SLiM data
        s_type, s_start, s_seq, s_end = slim
            
        # Try to find out if the SLiM is in one of the
        # disordered regions

        # Initialize the disordered flag to False
        is_disordered = False

        # For each disordered region
        for row in disorder_df.iterrows():
            
            # If the SLiM is included in the disordered region
            if start >= row["start"] and end <= row["end"]:

                # Turn the disordered flag on
                is_disordered = True

        # Save information about the SLiM
        results.append({"source" : "disprot",
                        "slim_type" : s_type,
                        "slim_start" : s_start,
                        "slim_sequence" : s_seq, 
                        "slim_end" : s_end,
                        "disordered_region" : str(is_disordered)})

    # Convert the results into a data frame and return it
    return pd.DataFrame(results)
