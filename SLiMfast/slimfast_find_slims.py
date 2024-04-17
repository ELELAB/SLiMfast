#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_find_slims.py
#
#    Find SLiMs of interest in availabl experimental structures,
#    AlphaFold structures, or homologous structures of the protein
#    of interest. Check whether the SLiMs are in disordered regions.
#    If requested, perform solvent accessibility analysis on residues
#    of interest in the SLiMs.
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


# Standard library
import json
import logging as log
import os
import sys
# Third-party packages
import pandas as pd
# SLiMfast
from .argparse_tools import get_arguments
from . import defaults, util


#######################################################################


def run():

    # Set the script's name
    name = "slimfast_find_slims"

    # Set the script's description
    description = \
        "\nFind SLiMs of interest in available experimental " \
        "structures, AlphaFold structures, or homologous " \
        "structures of the protein of interest. Check whether " \
        "the SLiMs are in disordered regions. If requested, " \
        "perform solvent accessibility analysis on residues " \
        "of interest in the SLiMs.\n"

    # Set the names of the relevant arguments to be parsed from the
    # command line
    args_names = \
        ["uniprot_id", "input_protein_file", "input_alphafold_file",
         "input_homologs_file", "input_alignments_file",
         "input_slims_file", "output_slims_protein_file",
         "output_slims_alphafold_file", "output_slims_homologs_file",
         "output_slims_disprot_file", "configrun_file",
         "log_file", "v", "vv", "log_to_stdout"]

    # Set the default values for the arguments
    default_values = \
        {"log_file" : defaults.LOG_FILES["slimfast_find_slims"]}

    #-----------------------------------------------------------------#

    # Add relevant arguments to the argument parser and return them
    args = get_arguments(name = name,
                         description = description,
                         args_names = args_names,
                         default_values = default_values)

    #-----------------------------------------------------------------#

    # Retrieve the arguments
    uniprot_id = args.uniprot_id
    input_protein_file = args.input_protein_file
    input_alphafold_file = args.input_alphafold_file
    input_homologs_file = args.input_homologs_file
    input_alignments_file = args.input_alignments_file
    input_slims_file = args.input_slims_file
    output_slims_protein_file = args.output_slims_protein_file
    output_slims_alphafold_file = args.output_slims_alphafold_file
    output_slims_homologs_file = args.output_slims_homologs_file
    output_slims_disprot_file = args.output_slims_disprot_file
    configrun_file = args.configrun_file
    log_file = args.log_file
    log_to_stdout = args.log_to_stdout
    v = args.verbose
    vv = args.verbose_debug

    #-----------------------------------------------------------------#

    # Get the default logging level
    log_level = defaults.LOG_LEVEL

    # If the user requested verbose logging
    if v:

        # The minimal logging level will be INFO
        log_level = log.INFO

    # If the user requested logging for debug purposes
    # (-vv overrides -v if both are provided)
    if vv:

        # The minimal logging level will be DEBUG
        log_level = log.DEBUG

    # Get the handlers
    handlers = \
        util.get_log_handlers(\
            # Whether to log to the standard output or not
            log_to_stdout = log_to_stdout,
            # The class of the file handler, if any
            log_file_handler_class = log.FileHandler,
            # The options to set the file handler
            log_file_handler_options = {"filename" : log_file,
                                        "mode" : "w"},
            # The level below which log messages are silenced
            log_level = log_level)

    # Set the logging configuration
    log.basicConfig(# The level below which log messages are silenced
                    level = log_level,
                    # The handlers
                    handlers = handlers,
                    # Other logging options
                    **defaults.LOG_OPTIONS)

    #-----------------------------------------------------------------#

    # Try to load the configuration for running
    try:
        
        config_run = \
            util.load_config(config_file = configrun_file,
                             config_type = "slimfast_find_slims")

    # If something went wrong 
    except Exception as e:
        
        # Set the error message
        errstr = \
            "It was not possible to load the configuration from " \
            f"'{configrun_file}'."

        # Log the exception
        log.exception(f"{errstr}. Error: {e}")

        # Re-raise the exception
        raise Exception(errstr) from e

    # Inform the user that the configuration was successfully loaded
    infostr = \
        "The configuration was successfully loaded from " \
        f"'{configrun_file}'."
    log.info(infostr)

    #-----------------------------------------------------------------#

    # Try to load the SLiM instances
    try:
        
        slims_to_find = util.load_slims(slims_file = input_slims_file)
    
    # If something went wrong
    except Exception as e:
        
        # Set the error message
        errstr = \
            "It was not possible to load the SLiM instances " \
            f"from '{input_slims_file}'."

        # Log the exception
        log.exception(f"{errstr}. Error: {e}")

        # Re-raise the exception
        raise Exception(errstr) from e

    # Inform the user that the SLiM instances were successfully
    # lodead
    infostr = \
        "The SLiM instances were successfully loaded " \
        f"from '{input_slims_file}'."
    log.info(infostr)

    #-----------------------------------------------------------------#

    # If the CSV file with the data frame with data about the
    # structures available in the Protein Data Bank for the protein of
    # interest was passed
    if input_protein_file is not None:

        # Get the field separator used in the data frame
        sep = defaults.OUTPUT_FILES["protein"]["options"]["sep"]
    
        # Try loading the data frame
        try:

            protein_df = pd.read_csv(input_protein_file,
                                     sep = sep)
        
        # If something went wrong
        except Exception as e:

            # Set the error message
            errstr = \
                "It was not possible to load the data frame " \
                "containing data about the structures available " \
                "in the Protein Data Bank for the protein of " \
                f"interest from '{input_protein_file}'."

            # Log the exception
            log.exception(f"{errstr}. Error: {e}")

            # Re-raise the exception
            raise Exception(errstr) from e

        # Inform the user that the data frame was successfully
        # loaded
        infostr = \
            "The data frame containing data about the structures " \
            "available in the Protein Data Bank for the protein " \
            "of interest was successfully loaded from " \
            f"'{input_protein_file}.'"
        log.info(infostr)

    # If the file was not passed
    else:

        # The corresponding data frame will be empty
        protein_df = pd.DataFrame()

    #-----------------------------------------------------------------#

    # If the CSV file with the data frame with data about the
    # structures available in the AlphaFold Protein Structure Database
    # for the protein of interest was passed
    if input_alphafold_file is not None:

        # Get the field separator used in the data frame
        sep = defaults.OUTPUT_FILES["alphafold"]["options"]["sep"]
        
        # Try loading the data frame
        try:
            
            alphafold_df = pd.read_csv(input_alphafold_file,
                                       sep = sep)
        
        # If something went wrong
        except Exception as e:

            # Set the error message
            errstr = \
                "It was not possible to load the data frame " \
                "containing data about the structures available in " \
                "the AlphaFold Protein Structure Database for the " \
                f"protein of interest from '{input_alphafold_file}'."

            # Log the exception
            log.exception(f"{errstr}. Error: {e}")

            # Re-raise the exception
            raise Exception(errstr) from e

        # Inform the user that the data frame was successfully
        # loaded
        infostr = \
            "The data frame containing data about the structures " \
            "available in the AlphaFold Protein Structure Database " \
            "for the protein of interest was successfully loaded " \
            f"from '{input_alphafold_file}'."
        log.info(infostr)

    # If the file was not passed
    else:

        # The corresponding data frame will be empty
        alphafold_df = pd.DataFrame()

    #-----------------------------------------------------------------#

    # If the CSV file with the data frame with data about the
    # structures of the homologs of the protein of interest available
    # in the Protein Data Bank was passed together with the JSON file
    # containing the corresponding alignments
    if input_homologs_file is not None \
    and input_alignments_file is not None:

        # Get the field separator used in the data frame
        sep = defaults.OUTPUT_FILES["homologs"]["options"]["sep"]

        # Try loading the data frame
        try:
            
            homologs_df = pd.read_csv(input_homologs_file,
                                      sep = sep)
        
        # If something went wrong
        except Exception as e:

            # Set the error message   
            errstr = \
                "It was not possible to load the data frame " \
                "containing data about the structures available " \
                "in the Protein Data Bank for the homologs of the " \
                f"protein of interest from '{input_homologs_file}'."

            # Log the exception
            log.exception(f"{errstr}. Error: {e}")

            # Re-raise the exception
            raise Exception(errstr) from e

        # Inform the user that the data frame was successfully
        # loaded
        infostr = \
            "The data frame containing data about the structures " \
            "available in the Protein Data Bank for the homologs " \
            "of the protein of interest was successfully loaded " \
            f"from '{input_homologs_file}'."
        log.info(infostr)

        # Try loading the alignments
        try:
            
            alignments = json.load(open(input_alignments_file, "r"))
        
        # If something went wrong
        except Exception as e:
            
            # Set the error message
            errstr = \
                "It was not possible to load the alignments from " \
                f"'{input_alignments_file}'."

            # Log the exception
            log.exception(f"{errstr}. Error: {e}")

            # Re-raise the exception
            raise Exception(errstr) from e

        # Inform the user that the alignments were successfully loaded
        infostr = \
            "The alignments were successfully loaded from " \
            f"'{input_alignments_file}'."
        log.info(infostr)

    # If the CSV file with the data frame with data about the
    # structures of the homologs of the protein of interest available
    # in the Protein Data Bank was not passed, but the JSON file
    # containing the corresponding alignments was
    elif input_homologs_file is None \
    and input_alignments_file is not None:

        # Set the error message
        errstr = \
            "If you pass '-ial', '--input-alignments-file', you " \
            "must also pass '-ih', '--input-homologs-file'."

        # Raise an exception
        raise Exception(errstr)

    # If the CSV file with the data frame with data about the
    # structures of the homologs of the protein of interest available
    # in the Protein Data Bank was passed, but the JSON file
    # containing the corresponding alignments was not
    elif input_homologs_file is not None \
    and input_alignments_file is None:

        # Set the error message
        errstr = \
            "If you pass '-ih', '--input-homologs-file', you " \
            "must also pass '-ial', '--input-alignments-file'."

        # Raise an exception
        raise Exception(errstr)

    # If neither the CSV file with the data frame with data about the
    # structures of the homologs of the protein of interest available
    # in the Protein Data Bank, nor the JSON file containing the
    # corresponding alignments were passed
    elif input_homologs_file is None \
    and input_alignments_file is None:

        # The corresponding data frame/dictionary will be empty
        homologs_df = pd.DataFrame()
        alignments = {}

    #-----------------------------------------------------------------#

    # Get the working directory where the program will be run
    wd = config_run["wd"]

    #-----------------------------------------------------------------#

    # Get the minimum length for a region to be considered
    # disordered
    min_len_disorder = config_run["structures"]["min_len_disorder"]

    # Get the name of the directory where to store the PDB files
    struct_dir = \
        os.path.join(wd, config_run["structures"]["struct_dir"])

    # Get the options for running DSSP
    dssp_params = config_run["dssp"]["run"]

    # Get the options for processing the DSSP output
    dssp_process_params = config_run["dssp"]["processing"]
    dssp_mapping = dssp_process_params["dssp_mapping"]
    disorder_symbols = dssp_process_params["disorder_symbols"]

    # Get the options to run NACCESS
    naccess_params = config_run["naccess"]["run"]

    # Create an empty list to store the structures and the SLiMs
    # found in them
    structures_and_slims = []

    #-----------------------------------------------------------------#

    # If there we have the structures available in the Protein Data
    # Bank for the protein of interest
    if not protein_df.empty:

        # Try to get the SLiMs found in each structure and the SLiMs
        # still to find
        try:
            
            structures_and_slims_prot, still_to_find = \
                util.find_slims(source = "structures",
                                structures_df = protein_df,
                                slims_to_find = slims_to_find,
                                struct_dir = struct_dir)

            # Add the structures and the SLiMs found in them
            # to the list
            structures_and_slims.extend(structures_and_slims_prot)
        
        # If something went wrong
        except Exception as e:

            # Set the error message
            errstr = \
                "It was not possible look for the SLiMs in the " \
                "structures available in the Protein Data Bank " \
                "for the protein of interest."

            # Re-raise the exception
            raise Exception(errstr) from e

        # Update the list of SLiMs still to be found
        slims_to_find = still_to_find

    #-----------------------------------------------------------------#

    # If there are still SLiMs to find and we have the structures
    # available in the AlphaFold Protein Structure Database for the
    # protein of interest
    if slims_to_find and not alphafold_df.empty:

        # Try to get the SLiMs found in each structure and the SLiMs
        # still to find
        try:
            
            structures_and_slims_alpha, still_to_find = \
                util.find_slims(source = "alphafold", 
                                structures_df = alphafold_df,
                                slims_to_find = slims_to_find,
                                struct_dir = struct_dir)

            # Add the structures and the SLiMs found in them
            # to the list
            structures_and_slims.extend(structures_and_slims_alpha)
        
        # If something went wrong
        except Exception as e:
            
            # Set the error message
            errstr = \
                "It was not possible to look for the SLiMs in the " \
                "structures available in the AlphaFold Protein " \
                "Structure Database for the protein of interest."

            # Re-raise the exception
            raise Exception(errstr) from e

        # Update the list of SLiMs still to be found
        slims_to_find = still_to_find

    #-----------------------------------------------------------------#

    # If there are still SLiMs to find and we have the structures
    # available in the Protein Data Bank for the homologs of the
    # protein of interst
    if slims_to_find and not homologs_df.empty:

        # Try to get the SLiMs found in each structure and the SLiMs
        # still to find
        try:
            
            structures_and_slims_hom, still_to_find = \
                util.find_slims(source = "homologs", 
                                structures_df = homologs_df,
                                slims_to_find = slims_to_find,
                                struct_dir = struct_dir,
                                als = alignments)

            # Add the structures and the SLiMs found in them
            # to the list
            structures_and_slims.extend(structures_and_slims_hom)
        
        # If something went wrong
        except Exception as e:
            
            # Set the error message
            errstr = \
                "It was not possible to look for the SLiMs in the " \
                "structures available in the Protein Data Bank for " \
                "the homologs of the protein of interest."

            # Re-raise the exception
            raise Exception(errstr) from e

        # Update the list of SLiMs still to be found
        slims_to_find = still_to_find

    #-----------------------------------------------------------------#

    # Try to get the results
    try:

        results_protein_df, \
        results_alphafold_df, \
        results_homologs_df = \
            util.get_results(\
                structures_and_slims = structures_and_slims,
                dssp_params = dssp_params,
                dssp_mapping = dssp_mapping,
                naccess_params = naccess_params,
                min_len_disorder = min_len_disorder,
                disorder_symbols = disorder_symbols,
                wd = wd)
    
    # If something went wrong
    except Exception as e:

        # Report it and exit
        errstr = \
            "It was not possible to create the data frames with " \
            "the SLiMs found in the available structures."

        # Re-raise the exception
        raise Exception(errstr) from e

    #-----------------------------------------------------------------#

    # If no SLiMs were found in the structures available in the Protein
    # Data Bank for the protein of interst
    if results_protein_df.empty:

        # If a data frame of structures was initially passed
        if not protein_df.empty:

            # Inform the user of it
            warnstr = \
                "No SLiMs were found in the structures available " \
                "in the Protein Data Bank for the protein of " \
                "interest."
            log.warning(warnstr)

    # Otherwise
    else:

        # Save the SLiMs found in the structures in the output
        # CSV file
        results_protein_df.to_csv(\
            os.path.join(wd, output_slims_protein_file),
            **defaults.OUTPUT_FILES["slims_protein"]["options"])

        # Inform the user that the data frame was successfully saved
        infostr = \
            "The data frame containing the SLiMs found in the " \
            "structures available in the Protein Data Bank for " \
            "the protein of interest was successfully saved in " \
            f"'{output_slims_protein_file}'."
        log.info(infostr)

    #-----------------------------------------------------------------#

    # If no SLiMs were found in the structures available in the
    # AlphaFold Protein Structure Database for the protein of interest
    if results_alphafold_df.empty:

        # If a data frame of structures was initially passed
        if not alphafold_df.empty:

            # Inform the user of it
            warnstr = \
                "No SLiMs were found in the structures available " \
                "in the AlphaFold Protein Structure Database for " \
                "the protein of interest."
            log.warning(warnstr)

    # Otherwise
    else:

        # Save the SLiMs found in the structures in the output
        # CSV file
        results_alphafold_df.to_csv(\
            os.path.join(wd, output_slims_alphafold_file),
            **defaults.OUTPUT_FILES["slims_alphafold"]["options"])

        # Inform the user that the data frame was successfully saved
        infostr = \
            "The data frame containing the SLiMs found in the " \
            "structures available in the AlphaFold Protein " \
            "Structure Database for the protein of interest was " \
            f"successfully saved in '{output_slims_alphafold_file}'."
        log.info(infostr)

    #-----------------------------------------------------------------#

    # If no SLiMs were found in the structures available in the
    # Protein Data Bank for the homologs of the protein of interest
    if results_homologs_df.empty:

        # If a data frame of structures was initially passed
        if not homologs_df.empty:

            # Inform the user of it
            warnstr = \
                "No SLiMs were found in the structures available " \
                "in the Protein Data Bank for the homologs of the " \
                "protein of interest."
            log.warning(warnstr)

    # Otherwise
    else:

        # Save the SLiMs found in the structures in the output
        # CSV file
        results_homologs_df.to_csv(\
            os.path.join(wd, output_slims_homologs_file),
            **defaults.OUTPUT_FILES["slims_homologs"]["options"])

        # Inform the user that the data frame was successfully saved
        infostr = \
            "The data frame containing the SLiMs found in the " \
            "structures available in the Protein Data Bank for " \
            "the homologs of the protein of interest was " \
            f"successfully saved in '{output_slims_homologs_file}'."
        log.info(infostr)

    #-----------------------------------------------------------------#

    # If there are still SLiMs that were not found
    if slims_to_find:

        # Try to get their DisProt annotations
        try:
            
            results_disprot_df = \
                util.get_disprot_annotations(\
                    uniprot_id = uniprot_id,
                    slims_to_find = slims_to_find)

        # If something went wrong
        except Exception as e:

            # Set the error message
            errstr = \
                "It was not possible to get the DisProt annotations " \
                "for the SLiMs of interest."

            # Re-raise the exception
            raise Exception(errstr) from e

        # If no DisProt annotations were found
        if results_disprot_df.empty:

            # Inform the user of it
            warnstr = \
                "No annotations were found in DisProt for the SLiMs " \
                "of interest."
            log.warning(warnstr)

        # Otherwise
        else:

            # Save the SLiMs annotated through DisProt in the output
            # CSV file
            results_disprot_df.to_csv(\
                os.path.join(wd, output_slims_disprot_file), 
                **defaults.OUTPUT_FILES["slims_disprot"]["options"])

            # Inform the user that the data frame was successfully
            # saved
            infostr = \
                "The data frame containing the SLiMs annotated " \
                "through DisProt was successfully saved in " \
                f"'{output_slims_disprot_file}'."
            log.info(infostr)


#######################################################################


def main():
    
    # Try to run
    try:

        run()

    # If something went wrong
    except Exception as e:

        # Report it and exit
        errstr = \
            "It was not possible to run " \
            f"'slimfast_find_slims'. Exception: {e}"
        log.exception(errstr)
        sys.exit(errstr)
