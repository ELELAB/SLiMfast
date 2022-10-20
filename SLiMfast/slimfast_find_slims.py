#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_find_slims.py
#
#    Find predicted SLiMs in a protein sequence in the available
#    PDB structures of the protein and/or in the available structures
#    of the protein's homologs, and check whether they belong to a
#    disordered region
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
import json
import logging as log
import os
import sys
# Third-party packages
import pandas as pd
# SLiMfast
from SLiMfast import util



# Get the module logger
logger = log.getLogger(__name__)



def get_args():


    #------------------------ Argument parser ------------------------#
    

    # Set the script description
    description = \
        "\nFind SLiMs in structures and check if they are in " \
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
                              log_to_stdout = True,
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
    log_to_stdout = args.log_to_stdout
    v = args.verbose
    vv = args.verbose_debug


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
    
    # Create the list containing the log handlers
    handlers = [log.FileHandler(filename = log_file, mode = "w")]

    # If the user selected the option to log also to the standard
    # output
    if log_to_stdout:

        # Add the StreamHandler to the list of handlers
        handlers.append(log.StreamHandler(sys.stdout))
    
    # Set the basic logging configuration
    log.basicConfig(level = level,
                    handlers = handlers)


    #------------------------- Configuration -------------------------#


    # Try to load the configuration file for running
    try:
        
        config_run = util.read_config_file(configrun_file)
    
    except Exception as e:
        
        # If something went wrong, report it and exit
        errstr = \
            f"Could not parse the configuration file " \
            f"{configrun_file}: {e}"
        logger.exception(errstr)
        sys.exit(errstr)

    # Inform the user of the successful parsing of the
    # configuration file
    infostr = \
        f"The configuration has been successfully parsed " \
        f"from {configrun_file}."
    logger.info(infostr)


    #-------------------------- Input files --------------------------#


    # Try to load the text file containing the SLiM instances
    try:
        
        slims_to_find = \
            util.read_slims_file(slims_file = input_slims_file)

        # Inform the user of the successful parsing of the SLiMs'
        # file
        infostr = \
            f"The SLiM instances have been successfully retrieved " \
            f"from {input_slims_file}."
        logger.info(infostr)
    
    # If something went wrong, report it and exit
    except Exception as e:
        
        errstr = f"Could not parse the file {input_slims_file}."
        logger.exception(errstr)
        sys.exit(errstr)

    # If the CSV file with the available structures of the
    # protein of interest has been passed
    if input_protein_file is not None:
        
        # Try to load the data frame containing the structures of
        # the protein of interest
        try:
            
            sep = config_run["structures"]["protein"]["out"]["sep"]
            protein_df = pd.read_csv(input_protein_file,
                                     sep = sep)

            # Inform the user of the successful parsing of the file
            # containing the available structures of the protein of
            # interest
            infostr = \
                f"Data about the available structures of the " \
                f"protein with UniProt ID {uniprot_id} have been " \
                f"successfully retrieved from {input_protein_file}."
            logger.info(infostr)
        
        # If something went wrong
        except Exception as e:

            # Report it and exit
            errstr = f"Could not parse the file {input_protein_file}."
            logger.exception(errstr)
            sys.exit(errstr)

    # If the CSV file with the available PDB structures of the
    # homologs has been passed
    if input_homologs_file is not None \
        and input_alignments_file is not None:

        # Try to load the data frame containing the structures of
        # the homologs of the protein of interest
        try:
            
            sep = config_run["structures"]["homolog"]["out"]["sep"]
            homologs_df = pd.read_csv(input_homologs_file,
                                      sep = sep)

            # Inform the user of the successful parsing of the file
            # containing the available structures of the homologs of
            # the protein of interest
            infostr = \
                f"Data about the available structures of the " \
                f"homologs of the protein with UniProt ID {uniprot_id} " \
                f"have been successfully retrieved from " \
                f"{input_homologs_file}."
            logger.info(infostr)
        
        # If something went wrong
        except Exception as e:

            # Report it and exit     
            errstr = f"Could not parse the file {input_homologs_file}."
            logger.exception(errstr)
            sys.exit(errstr)

        # Try to load the alignments
        try:
            
            alignments = json.load(open(input_alignments_file, "r"))

            # Inform the user of the successful parsing of the file
            # containing the alignments of the homologs to the
            # protein of interest
            infostr = \
                f"Data about the alignments of the homologs to the " \
                f"protein with UniProt ID {uniprot_id} have been " \
                f"successfully retrieved from {input_alignments_file}."
            logger.info(infostr)
        
        # If something went wrong
        except Exception as e:
            
            # Report it and exit
            errstr = \
                f"Could not parse the file {input_alignments_file}."
            logger.exception(errstr)
            sys.exit(errstr)

    # If the CSV of the homologs has been passed but not the
    # alignments' JSON file
    elif input_homologs_file is None \
        and input_alignments_file is not None:

            # Report it and exit
            errstr = \
                f"If you pass the JSON file containing the " \
                f"alignments of the homologs of the protein with " \
                f"UniProt ID {uniprot_id} with the protein itself, " \
                f"you must also pass the CSV file containing the " \
                f"data about the available structures of the homologs."
            logger.exception(errstr)
            sys.exit(errstr)

    # If the alignments' JSON file has been passed but not the
    # homologs' CSV file
    elif input_homologs_file is not None \
        and input_alignments_file is None:

            # Report it and exit
            errstr = \
                f"If you pass the CSV file containing the data " \
                f"about the available structures of the homologs " \
                f"of the protein with UniProt ID {uniprot_id}, you " \
                f"must also pass the JSON file containing the " \
                f"alignments of the homologs to the protein."
            logger.exception(errstr)
            sys.exit(errstr)

    # Return data frames, alignments, SLiM instances, output files
    # and configuration
    return uniprot_id, protein_df, homologs_df, alignments, \
           slims_to_find, output_slims_protein_file, \
           output_slims_homologs_file, output_slims_disprot_file, \
           config_run


def run(uniprot_id,
        protein_df,
        homologs_df,
        alignments,
        slims_to_find,
        output_slims_protein_file,
        output_slims_homologs_file,
        output_slims_disprot_file,
        config_run):
 

    # Get the working directory where the pipeline will be run
    wd = util.get_wd(config_run)

    # Get the minimum length for a region to be considered
    # disordered and the name of the directory where to store
    # the PDB files
    min_len_disorder = config_run["structures"]["min_len_disorder"]
    struct_dir = \
        os.path.join(wd, config_run["structures"]["struct_dir"])

    # Get the options for running DSSP
    dssp_params = config_run["dssp"]["run"]

    # Get the options for processing the DSSP output
    dssp_process_params = config_run["dssp"]["processing"]
    dssp_mapping = dssp_process_params["dssp_mapping"]
    disorder_symbols = dssp_process_params["disorder_symbols"]

    # Get the parameters for the output data frames
    out_params = config_run["slims"]["out"]

    structures_and_slims = []


    #------------------- Find SLiMs in structures --------------------#


    # If there are structures available for the protein of interest
    if protein_df is not None:

        # Report the SLiMs found in each PDB structure and the SLiMs
        # still to find
        try:
            
            structures_and_slims_prot, still_to_find = \
                util.find_slims(source = "protein",
                                structures_df = protein_df,
                                slims_to_find = slims_to_find,
                                struct_dir = struct_dir)

            structures_and_slims.extend(structures_and_slims_prot)
        
        # If something went wrong
        except Exception as e:

            # Report it and exit
            errstr = \
                f"Could not search for the SLiMs in the structures " \
                f"of the protein with UniProt ID {uniprot_id}: {e}"
            logger.exception(errstr)
            sys.exit(errstr)

        # Update the list of SLiMs still to be found
        slims_to_find = still_to_find


    #-------------------- Find SLiMs in homologs ---------------------#


    # If there are no structures available for the given protein
    # or some SLiMs were not found in the available structures
    if homologs_df is not None and slims_to_find:

        # Report the SLiMs found in each structure and the SLiMs
        # still to find
        try:
            
            structures_and_slims_hom, still_to_find = \
                util.find_slims(source = "homologs", 
                                structures_df = homologs_df,
                                slims_to_find = slims_to_find,
                                struct_dir = struct_dir,
                                als = alignments)

            structures_and_slims.extend(structures_and_slims_hom)
        
        # If something went wrong
        except Exception as e:
            
            # Report it and exit
            errstr = \
                f"Could not search for the SLiMs in the structures " \
                f"of the homologs of the protein with UniProt " \
                f"ID {uniprot_id}: {e}"
            logger.exception(errstr)
            sys.exit(errstr)

        # Update the list of SLiMs still to be found
        slims_to_find = still_to_find

    # Try to get the results
    try:

        results_protein_df, results_homologs_df = \
            util.get_results(\
                structures_and_slims = structures_and_slims,
                dssp_params = dssp_params,
                dssp_mapping = dssp_mapping,
                min_len_disorder = min_len_disorder,
                disorder_symbols = disorder_symbols,
                wd = wd)
    
    # If something went wrong
    except Exception as e:

        # Report it and exit
        errstr = \
            f"Could not gather the results for the SLiMs found. " \
            f"Exception: {e}"
        logger.exception(errstr)
        sys.exit(errstr)

    # If no SLiMs were found in the structures of the protein of
    # interest, and the original data frame had available
    # structures for the protein of interest
    if results_protein_df.empty and (not protein_df.empty):

        # Inform the user of it
        warnstr = \
            f"No SLiMs found in the structures of the protein with " \
            f"UniProt ID {uniprot_id}."
        logger.warning(warnstr)

    # Otherwise
    else:

        # Write out the SLiMs found in the structures of the protein
        # of interest
        results_protein_df.to_csv(\
            os.path.join(wd, output_slims_protein_file), **out_params)

        # Inform the user of the successful writing of the output file
        infostr = \
            f"The SLiMs found in the structures of the protein with " \
            f"UniProt ID {uniprot_id} have been successfully " \
            f"written to {output_slims_protein_file}."
        logger.info(infostr)
    
    # If no SLiMs were found in the structures of the homologs of the
    # protein of interest, and the original data frame had available
    # structures for the homologs, and no SLiMs were found in the
    # available structures of the protein
    if results_homologs_df.empty and (not homologs_df.empty) \
    and results_protein_df.empty:

        # Inform the user of it
        warnstr = \
            f"No SLiMs found in the structures of the homologs of " \
            f"the protein with UniProt ID {uniprot_id}."
        logger.warning(warnstr)

    # Otherwise
    else:

        # Write out the SLiMs found in the structures of the homologs
        # of the protein of interest
        results_homologs_df.to_csv(\
            os.path.join(wd, output_slims_homologs_file), **out_params)

        # Inform the user of the successful writing of the output file
        infostr = \
            f"The SLiMs found in the structures of the homologs of " \
            f"the protein with UniProt ID {uniprot_id} have been " \
            f"successfully written to {output_slims_protein_file}."
        logger.info(infostr)


    #-------------------- Get DisProt annotations --------------------#


    # If there are still SLiMs that were not annotated
    if slims_to_find:

        # Try to find their DisProt annotations
        try:
            results_disprot_df = \
                util.get_disprot_annotations(\
                    uniprot_id = uniprot_id,
                    slims_to_find = slims_to_find)

        # If something went wrong
        except Exception as e:

            # Report it and exit
            errstr = \
                f"Could not get the DisProt annotations for the " \
                f"SLiMs of interest in the protein with UniProt ID " \
                f"{uniprot_id}. Exception: {e}"
            logger.exception(errstr)
            sys.exit(errstr)

        # If no DisProt annotations were found
        if results_disprot_df.empty:

            # Inform the user of it
            warnstr = \
                f"No annotation found in DisProt for the SLiMs " \
                f"of interest in the protein with UniProt ID " \
                f"{uniprot_id}."
            logger.warning(warnstr)

        # Otherwise
        else:

            # Write out the SLiMs annotated through DisProt
            results_disprot_df.to_csv(\
                os.path.join(wd, output_slims_disprot_file), 
                             **out_params)

            # Inform the user of the successful writing of the output file
            infostr = \
                f"The SLiMs annotated through DisProt for the " \
                f"protein with UniProt ID {uniprot_id} have been " \
                f"successfully written to {output_slims_protein_file}."
            logger.info(infostr)



def main():
    args = get_args()
    run(*args)