#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_find_slims.py
#
#    Find predicted SLiMs in a protein sequence in the available
#    PDB structures of the protein and/or in the available structures
#    of the protein's homologs, and check whether they belong to a
#    disordered region
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
import copy
import logging as log
import os
import sys
# Third-party packages
import yaml
# SLiMfast
from SLiMfast import util



def get_args():


    #------------------------ Argument parser ------------------------#
    

    # Set the script description
    description = \
        "\nFind SLiMs in PDB structures and check if they are in " \
        "disordered regions.\n"

    # Set the default name for the log file
    log_file_default = "slimfast_find_slims.log"

    # Add relevant arguments to the argument parser and return them
    args = util.get_arguments(description = description,
                              uniprot_id = True,
                              input_protein_file = True,
                              input_homologs_file = True,
                              input_alignments_file = True,
                              input_slims_file = True,
                              output_slims_protein_file = True,
                              output_slims_homologs_file = True,
                              output_slims_disprot_file = True,
                              configrun_file = True,
                              log_file = True,
                              v = True,
                              vv = True,
                              log_file_default = log_file_default)

    # Retrieve the arguments
    uniprot_id = args.uniprot_id
    input_protein_file = args.input_protein_file
    input_homologs_file = args.input_homologs_file
    input_alignments_file = args.input_alignments_file
    input_slims_file = args.input_slims_file
    output_slims_protein_file = args.output_slims_protein_file
    output_slims_homologs_file = args.output_slims_homologs_file
    output_slims_disprot_file = args.output_slims_disprot_file
    configrun_file = args.configrun_file
    log_file = args.log_file
    v = args.v
    vv = args.vv


    #---------------------------- Logging ----------------------------#


    # Default logging level is WARNING
    level = log.WARNING
    # If verbose logging required
    if v:
        # Logging level will be INFO
        level = log.INFO
    # If debug logging required
    if vv:
        # Logging level will be DEBUG
        level = log.DEBUG
    
    # Basic logging configuration
    log.basicConfig(level = level, 
                    filename = log_file, 
                    filemode = "w")


    #------------------------- Configuration -------------------------#


    # Try to load the configuration file for running
    try:
        config_run = yaml.safe_load(open(configrun_file, "r"))
    except Exception as e:
        # If something went wrong, report it and exit
        errstr = f"Could not parse the run configuration file " \
                 f"{configrun_file}: {e}"
        log.error(errstr)
        sys.exit(errstr)

    # Try to load the text file containing the SLiM instances
    try:
        slim_instances = \
            util.read_slim_instances(slim_file = input_slims_file)
    # If something went wrong, report it and exit
    except Exception as e:
        errstr = f"Could not parse the file {input_slims_file}."
        log.error(errstr)
        sys.exit(errstr)


    # Return dataframes, alignments, SLiM instances, output files
    # and configuration
    return protein_df, homologs_df, alignments, slims_instances, \
           output_slims_protein_file, output_slims_homologs_file, \
           config_run


def run(protein_df,
        homologs_df,
        alignments,
        slims_instances,
        output_slims_protein_file,
        output_slims_homologs_file,
        config_run):
 

    # Get the working directory
    wd = config_run["wd"]

    # Get whether to find the core or extended SLiM in the
    # structures
    find_in_pdb = struct_params["structures"]["find_in_pdb"]

    # Get the options for running DSSP
    dssp_params = config_run["dssp"]["run"]

    # Get the options for processing the DSSP output
    dssp_process_params = config_run["dssp"]["processing"]
    dssp_mapping = dssp_process_params["dssp_mapping"]
    min_len_disorder = dssp_process_params["min_len_disorder"]
    disorder_symbols = dssp_process_params["disorder_symbols"]
    
    # Make a copy of the original list of SLiMs to find (will
    # be modified)
    slims_to_find = copy.deepcopy(slims_instances)


    #------------------- Find SLiMs in structures --------------------#


    # If there are structures available for the protein of interest
    if protein_df:

        # Report the SLiMs found in each PDB structure and the SLiMs
        # still to find
        try:
            structures_and_slims, still_to_find = \
                util.find_slims(source = "protein",
                                structures_df = protein_df,
                                slims_to_find = slims_to_find,
                                find_in_pdb = find_in_pdb,
                                wd = wd)
        
        # If something went wrong, report it and exit
        except Exception as e:
            errstr = f"Could not find SLiMs in the PDB structures " \
                     f"of the protein of interest: {e}"
            log.error(errstr)
            sys.exit(errstr)

        # Update the list of SLiMs still to be found
        slims_to_find = still_to_find


    #-------------------- Find SLiMs in homologs ---------------------#


    # If there are no structures available for the given protein
    # or some SLiMs were not found in the available structures
    if homologs_df and (not protein_df or \
                        (protein_df and slims_to_find)):

        # Report the SLiMs found in each PDB structure and the SLiMs
        # still to find
        try:
            structures_and_slims, still_to_find = \
                util.find_slims(source = "homologs", 
                                structures_df = homologs_df,
                                als = alignments,
                                slims_to_find = slims_to_find,
                                find_in_pdb = find_in_pdb,
                                wd = wd)
        
        # If something went wrong, report it and exit
        except Exception as e:
            errstr = f"Could not find SLiMs in the PDB structures " \
                     f"of the homologs of the protein of interest: " \
                     f"{e}"
            log.error(errstr)
            sys.exit(errstr)

        # Update the list of SLiMs still to be found
        slims_to_find = still_to_find

    # Try to get the results
    try:
        results_protein_df, results_homologs_df = \
            util.get_results(\
                structures_and_slims = structures_and_slims,
                wd = wd)
    
    # If something went wrong, report it and exit
    except Exception as e:
        errstr = \
            f"Could not get the results for the SLiMs. Exception: {e}"
        logger.error(errstr)
        sys.exit(errstr)


    # Write out the CSV files containing the dataframes
    results_protein_df.to_csv(output_slims_protein_file,
                              **out_params)
    results_homologs_df.to_csv(output_slims_homologs_file, 
                               **out_params)


    #-------------------- Get DisProt annotations --------------------#


    # If there are still SLiMs that were not annotated
    if slims_to_find:

        # Try to find their DisProt annotations
        try:
            results_disprot_df = \
                util.get_disprot_annotations(\
                    uniprot_id = uniprot_id,
                    slims_to_find = slims_to_find,
                    wd = wd)

        # If something went wrong, report it and exit
        except Exception as e:
            errstr = \
                f"Could not get the DisProt annotations for the " \
                f"SLiMs. Exception: {e}"
            logger.error(errstr)
            sys.exit(errstr)


    # Write out the CSV file containing the dataframe
    results_disprot_df.to_csv(output_slims_disprot_file,
                              **out_params)



def main():
    args = get_args()
    run(*args)