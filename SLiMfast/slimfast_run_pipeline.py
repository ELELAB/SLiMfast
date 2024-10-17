#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_run_pipeline.py
#
#    Run the entire SLiMfast pipeline (or parts of it) for a series
#    of UniProt IDs in parallel.
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
import logging as log
import os
import sys
# Third-party packages
import dask
# SLiMfast
from .argparse_tools import get_arguments
from .run_tools import run_pipeline_step
from . import defaults, util


#######################################################################



def run():

    # Set the script's name
    name = "slimfast_run_pipeline"

    # Set the script's description
    description = \
        "\nRun the SLiMfast pipeline (or parts of it) in parallel " \
        "for multiple UniProt IDs.\n"

    # Set the names of the relevant arguments to be parsed from the
    # command line
    args_names = \
        ["input_ids_slims_file", "do_get_structures",
         "do_get_alphafold", "do_get_homologs", "do_find_slims",
         "output_protein_file", "output_alphafold_file",
         "output_homologs_file", "output_alignments_file",
         "output_slims_protein_file", "output_slims_alphafold_file",
         "output_slims_homologs_file", "output_slims_disprot_file",
         "configrun_file", "log_file", "v", "vv", "log_to_stdout",
         "nproc"]

    # Set the default values for the arguments
    default_values = \
        {"log_file" : defaults.LOG_FILES["slimfast_run_pipeline"]}

    #-----------------------------------------------------------------#

    # Add relevant arguments to the argument parser and return them
    args = get_arguments(name = name,
                         description = description,
                         args_names = args_names,
                         default_values = default_values)

    #-----------------------------------------------------------------#


    # Retrieve the arguments
    input_ids_slims_file = args.input_ids_slims_file
    do_get_structures = args.do_get_structures
    do_get_alphafold = args.do_get_alphafold
    do_get_homologs = args.do_get_homologs
    do_find_slims = args.do_find_slims
    output_protein_file = args.output_protein_file
    output_alphafold_file = args.output_alphafold_file
    output_homologs_file = args.output_homologs_file
    output_alignments_file = args.output_alignments_file
    output_slims_protein_file = args.output_slims_protein_file
    output_slims_alphafold_file = args.output_slims_alphafold_file
    output_slims_homologs_file = args.output_slims_homologs_file
    output_slims_disprot_file = args.output_slims_disprot_file
    configrun_file = util.get_abspath(args.configrun_file)
    log_file = args.log_file
    log_to_stdout = args.log_to_stdout
    v = args.verbose
    vv = args.verbose_debug
    nproc = args.nproc

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
                             config_type = "slimfast_run_pipeline")

    # If something went wrong 
    except Exception as e:
        
        # Set the error message
        errstr = \
            "It was not possible to load the configuration from " \
            f"'{configrun_file}'."

        # Re-raise the exception
        raise Exception(errstr) from e

    # Inform the user that the configuration was successfully loaded
    infostr = \
        "The configuration was successfully loaded from " \
        f"'{configrun_file}'."
    log.info(infostr)

    #-----------------------------------------------------------------#

    # Get the working directory where the pipeline will be run
    wd = config_run["wd"]

    #-----------------------------------------------------------------#

    # Get the logging configuration for Dask
    dask_logging_config = \
        util.get_dask_logging_config(log_to_stdout = log_to_stdout,
                                     log_file = log_file)

    #-----------------------------------------------------------------#
    
    # Set the configuration for Dask-specific logging
    dask.config.set({"distributed.logging" : dask_logging_config})

    # Import Dask-relates stuff here; otherwise, the logging does
    # not get configured properly
    from distributed import Client, LocalCluster

    # Create the local cluster
    cluster = LocalCluster(n_workers = nproc,
                           silence_logs = log.CRITICAL,
                           processes = True,
                           threads_per_worker = 1)
    
    # Open the client from the cluster
    client = Client(cluster)

    #-----------------------------------------------------------------#

    # Try to load the UniProt IDs and associated SLiMs to find
    try:

        ids_slims = util.load_ids_slims(input_ids_slims_file)

    # If something went wrong 
    except Exception as e:
        
        # Set the error message
        errstr = \
            "It was not possible to load the UniProt IDs and " \
            f"the associated SLiMs from '{input_ids_slims_file}'."

        # Re-raise the exception
        raise Exception(errstr) from e

    # Inform the that the UniProt IDs and associated SLiMs were
    # successfully loaded
    infostr = \
        "The UniProt IDs and associated SLiMs were successfully " \
        f"loaded from '{input_ids_slims_file}'."
    log.info(infostr)

    #-----------------------------------------------------------------#

    # Inform the user that the run is starting
    log.info("Running...")

    #-----------------------------------------------------------------#

    # Create a list of futures to be gathered in the end
    futures = []

    #-----------------------------------------------------------------#

    # For each UniProt ID and associated SLiMs
    for uniprot_id, slims_to_find in ids_slims.items():

        # Set the path to a directory with the UniProt ID inside the
        # working directory, where the results for the protein with
        # the current UniProt will be stored
        curr_wd = os.path.join(wd, uniprot_id)

        # Create the directory, if it does not exist
        os.makedirs(curr_wd,
                    exist_ok = True)

        #-------------------------------------------------------------#

        # Initialize the output CSV file for 'slimfast_get_structures'
        # to None
        protein_file = None

        # Initialize the output CSV file for 'slimfast_get_alphafold'
        # to None
        alphafold_file = None

        # Initialize the output CSV file for 'slimfast_get_homologs'
        # to None
        homologs_file = None

        # Initialize the output JSON file for 'slimfast_get_homologs'
        # to None
        alignments_file = None

        #-------------------------------------------------------------#

        # If the user requested the search for structures of the
        # protein of interest in the Protein Data Bank
        if do_get_structures:

            # Set the name of the output CSV file for 
            # 'slimfast_get_structures'
            protein_file = f"{uniprot_id}_{output_protein_file}"

        #-------------------------------------------------------------#

        # If the user requested the search for structures of the
        # protein of interest in the AlphaFold Protein Structure
        # Database
        if do_get_alphafold:

            # Set the name of the output CSV file for
            # 'slimfast_get_alphafold'
            alphafold_file = f"{uniprot_id}_{output_alphafold_file}"

        #-------------------------------------------------------------#

        # If the user requested the search for structures of the
        # homologs of the protein of interest in the Protein Data
        # Bank
        if do_get_homologs:

            # Set the name of the output CSV file for
            # 'slimfast_get_homologs'
            homologs_file = f"{uniprot_id}_{output_homologs_file}"

            # Set the name of the output JSON file for
            # 'slimfast_get_homologs'
            alignments_file = f"{uniprot_id}_{output_alignments_file}"

        #-------------------------------------------------------------#

        # If the user requested the SLiMs' search
        if do_find_slims:

            # Set the name of the input file containing the SLiMs
            input_slims_file_name = \
                defaults.INPUT_SLIMS_FILE["name"]

            # Set the options to write the input file containing the
            # SLiMs
            input_slims_file_options = \
                defaults.INPUT_SLIMS_FILE["options"]

            # Set the path to the SLiMs' input file
            slims_input_file = \
                os.path.join(curr_wd, 
                             f"{uniprot_id}_{input_slims_file_name}")

            # Write the SLiMs' input file
            slims_to_find.to_csv(slims_input_file,
                                 **input_slims_file_options)

            # Set the names of the output CSV files
            slims_protein_file = \
                f"{uniprot_id}_{output_slims_protein_file}"
            slims_alphafold_file = \
                f"{uniprot_id}_{output_slims_alphafold_file}"
            slims_homologs_file = \
                f"{uniprot_id}_{output_slims_homologs_file}"
            slims_disprot_file = \
                f"{uniprot_id}_{output_slims_disprot_file}"

            # Create a dictionary of temporary futures that will serve
            # as dependencies for 'slimfast_find_slims', if needed
            tmp_futures = {}

        #-------------------------------------------------------------#

        # If the user requested the search for structures of the
        # protein of interest in the Protein Data Bank
        if do_get_structures:

            # Get the name of the log file for this process
            process_log_file = \
                defaults.LOG_FILES["slimfast_get_structures"]

            # Launch the process
            output_files = \
                client.submit(\
                    run_pipeline_step,
                    step = "get_structures",
                    wd = curr_wd,
                    uniprot_id = uniprot_id,
                    output_files = [protein_file],
                    configrun_file = configrun_file,
                    log_file = process_log_file,
                    v = v,
                    vv = vv)

            # FOR FUTURE REFERENCE: the if/else block here is needed
            # because appending to the 'futures' list in any case
            # (i.e., whether 'do_find_slims' is True or False)
            # results in Dask throwing an error when 'do_find_slims'
            # is True. It is possible this happens because dask
            # tries to gather futures which have already completed
            # (since, for each UniProt ID, slimfast_get_structures
            # needs to complete before slimfast_find_slims is run,
            # if both are requested, meaning the future representing
            # the slimfast_get_structures process has already completed
            # by the end of the program).

            # If 'slimfast_find_slims' will be run
            if do_find_slims:

                # Add the future to the dictionary storing the
                # futures 'slimfast_find_slims' depends on
                tmp_futures["gs"] = output_files

            # Otherwise
            else:
                
                # Add the future to the list of futures that will
                # be gathered in the end
                futures.append(output_files)

        #-------------------------------------------------------------#

        # If the user requested the search for structures of the
        # protein of interest in the AlphaFold Protein Structure
        # Database
        if do_get_alphafold:

            # Get the name of the log file for this process
            process_log_file = \
                defaults.LOG_FILES["slimfast_get_alphafold"]

            # Launch the process
            output_files = \
                client.submit(\
                    run_pipeline_step,
                    step = "get_alphafold",
                    wd = curr_wd,
                    uniprot_id = uniprot_id,
                    output_files = [alphafold_file],
                    configrun_file = configrun_file,
                    log_file = process_log_file,
                    v = v,
                    vv = vv)

            # FOR FUTURE REFERENCE: the if/else block here is needed
            # because appending to the 'futures' list in any case
            # (i.e., whether 'do_find_slims' is True or False)
            # results in Dask throwing an error when 'do_find_slims'
            # is True. It is possible this happens because dask
            # tries to gather futures which have already completed
            # (since, for each UniProt ID, slimfast_get_alphafold
            # needs to complete before slimfast_find_slims is run,
            # if both are requested, meaning the future representing
            # the slimfast_get_alphafold process has already completed
            # by the end of the program).

            # If slimfast_find_slims will be run
            if do_find_slims:

                # Add the future to the dictionary storing the
                # futures slimfast_find_slims depends on
                tmp_futures["ga"] = output_files

            # Otherwise
            else:
                
                # Add the future to the list of futures that will
                # be gathered in the end
                futures.append(output_files)

        #-------------------------------------------------------------#

        # If the user requested the search for structures of the
        # homologs of the protein of interest in the Protein Data
        # Bank
        if do_get_homologs:

            # Get the name of the log file for this process
            process_log_file = \
                defaults.LOG_FILES["slimfast_get_homologs"]

            # Launch the process
            output_files = \
                client.submit(\
                    run_pipeline_step,
                    step = "get_homologs",
                    wd = curr_wd,
                    uniprot_id = uniprot_id,
                    output_files = [homologs_file, alignments_file],
                    configrun_file = configrun_file,
                    log_file = process_log_file,
                    v = v,
                    vv = vv)

            # FOR FUTURE REFERENCE: the if/else block here is needed
            # because appending to the 'futures' list in any case
            # (i.e., whether 'do_find_slims' is True or False)
            # results in Dask throwing an error when 'do_find_slims'
            # is True. It is possible this happens because dask
            # tries to gather futures which have already completed
            # (since, for each UniProt ID, slimfast_get_homologs
            # needs to complete before slimfast_find_slims is run,
            # if both are requested, meaning the future representing
            # the slimfast_get_homologs process has already completed
            # by the end of the program).

            # If slimfast_find_slims will be run
            if do_find_slims:

                # Add the future to the dictionary storing the
                # futures slimfast_find_slims depends on
                tmp_futures["gh"] = output_files

            # Otherwise
            else:
                
                # Add the futures to the list of futures that will
                # be gathered in the end
                futures.append(output_files)

        #-------------------------------------------------------------#

        # If the user requested the SLiMs' search
        if do_find_slims:

            # Get the input files
            input_files = \
                client.submit(\
                     util.get_values_from_dict,
                        d = tmp_futures,
                        keys = ["gs", "ga", "gh"],
                        expected_num_values = [1, 1, 2],
                        insert_first = [slims_input_file])

            # Set the output files for 'slimfast_find_slims'
            output_files = \
                [slims_protein_file, slims_alphafold_file,
                 slims_homologs_file, slims_disprot_file]

            # Get the name of the log file for this process
            process_log_file = \
                defaults.LOG_FILES["slimfast_find_slims"]

            # Launch the process
            output_files = \
                client.submit(\
                    run_pipeline_step,
                    step = "find_slims",
                    wd = curr_wd,
                    uniprot_id = uniprot_id,
                    input_files = input_files,
                    output_files = output_files,
                    configrun_file = configrun_file,
                    log_file = process_log_file,
                    v = v,
                    vv = vv)

            # Add the future to the list of futures that will
            # be gathered in the end
            futures.append(output_files)
    
    #-----------------------------------------------------------------#

    # Gather the remaining futures
    client.gather(futures)


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
            f"'slimfast_run_pipeline'. Exception: {e}"
        log.exception(errstr)
        sys.exit(errstr)
