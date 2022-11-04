#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_run_pipeline.py
#
#    Run the entire SLiMfast pipeline (or parts of it) for a series
#    of UniProt IDs in parallel.
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
import os
import sys
# Third-party packages
import dask
from distributed import (
    Client,
    LocalCluster
    )
# SLiMfast
from .argparse_tools import get_arguments
from .defaults import (
    INPUT_SLIMS_FILE,
    INPUT_SLIMS_FILE_OPTS,
    LOG_GET_STRUCTURES_FILE,
    LOG_GET_ALPHAFOLD_FILE,
    LOG_GET_HOMOLOGS_FILE,
    LOG_FIND_SLIMS_FILE,
    LOG_RUN_PIPELINE_FILE,
    )
from .run_tools import run_pipeline_step
from . import util



# Get the module logger
logger = log.getLogger(__name__)



def run():



    #------------------------ Argument parser ------------------------#
    


    # Set the script description
    description = \
        "\nRun the SLiMfast pipeline (or parts of it) in parallel " \
        "for multiple UniProt IDs.\n"

    # Add relevant arguments to the argument parser and return them
    args = get_arguments(\
                description = description,
                input_ids_slims_file = True,
                do_get_structures = True,
                do_get_alphafold = True,
                do_get_homologs = True,
                do_find_slims = True,
                output_protein_file = True,
                output_alphafold_file = True,
                output_homologs_file = True,
                output_alignments_file = True,
                output_slims_protein_file = True,
                output_slims_alphafold_file = True,
                output_slims_homologs_file = True,
                output_slims_disprot_file = True,
                configrun_file = True,
                log_file = True,
                log_to_stdout = True,
                v = True,
                vv = True,
                log_file_default = LOG_RUN_PIPELINE_FILE,
                nproc = True)

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

    # Get the working directory where the pipeline should be run
    wd = util.get_abspath(config_run["wd"])



    #---------------------- Cluster and client -----------------------#



    # Create the local cluster
    cluster = LocalCluster(n_workers = nproc,
                           silence_logs = log.CRITICAL,
                           processes = True,
                           threads_per_worker = 1)
    
    # Open the client from the cluster
    client = Client(cluster)



    #------------------ Input UniProt IDs and SLiMs ------------------#



    # Try to parse the CSV file containing the UniProt IDs and
    # corresponding SLiMs to find
    try:

        ids_slims = \
            util.read_ids_slims_file(input_ids_slims_file)

    # If something went wrong 
    except Exception as e:
        
        # Report it and exit
        errstr = \
            f"Could not parse the CSV file containing the UniProt " \
            f"IDs and corresponding SLiMs {input_ids_slims_file}: " \
            f" {e}"
        logger.exception(errstr)
        sys.exit(errstr)

    # Inform the user about the successful parsing of the file
    infostr = \
        f"UniProt IDs and corresponding SLiMs retrieved " \
        f"successfully from {input_ids_slims_file}."
    logger.info(infostr)



    #------------------------------ Run ------------------------------#



    # Inform the user that the run is starting
    logger.info("Running...")

    # Create a list of futures to be gathered in the end
    futures = []


    # For each UniProt ID and associated SLiMs
    for uniprot_id, slims_to_find in ids_slims.items():


        # Create a directory with the UniProt ID inside the
        # working directory, where the results for the current
        # protein will be stored
        curr_wd = os.path.join(wd, uniprot_id)
        os.makedirs(curr_wd, exist_ok = True)

        # Initialize the input files for slimfast_get_structures,
        # slimfast_get_alphafold, slimfast_get_homologs to None
        protein_file = None
        alphafold_file = None
        homologs_file = None
        alignments_file = None


        # If the user requested the search for structures of the
        # protein of interest in the Protein Data Bank
        if do_get_structures:

            # Set the name of the output/input file for
            # slimfast_get_structures/slimfast_find_slims
            protein_file = f"{uniprot_id}_{output_protein_file}"


        # If the user requested the search for structures of the
        # protein of interest in the AlphaFold Protein Structure
        # Database
        if do_get_alphafold:

            # Set the name of the output/input file for
            # slimfast_get_alphafold/slimfast_find_slims
            alphafold_file = f"{uniprot_id}_{output_alphafold_file}"


        # If the user requested the search for structures of the
        # homologs of the protein of interest in the Protein Data
        # Bank
        if do_get_homologs:

            # Set the name of the output/input files for
            # slimfast_get_homologs/slimfast_find_slims
            homologs_file = f"{uniprot_id}_{output_homologs_file}"
            alignments_file = f"{uniprot_id}_{output_alignments_file}"


        # If the user requested the SLiMs' search
        if do_find_slims:

            # Set the path to the SLiMs' input file
            slims_input_file = \
                os.path.join(curr_wd, 
                             f"{uniprot_id}_{INPUT_SLIMS_FILE}")

            # Write the SLiMs' input file
            slims_to_find.to_csv(slims_input_file,
                                 **INPUT_SLIMS_FILE_OPTS)

            # Set the names of the output files
            slims_protein_file = \
                f"{uniprot_id}_{output_slims_protein_file}"
            slims_alphafold_file = \
                f"{uniprot_id}_{output_slims_alphafold_file}"
            slims_homologs_file = \
                f"{uniprot_id}_{output_slims_homologs_file}"
            slims_disprot_file = \
                f"{uniprot_id}_{output_slims_disprot_file}"

            # Create a dictionary of temporary futures that will serve
            # as dependencies for slimfast_find_slims, if needed
            tmp_futures = {}


        # If the user requested the search for structures of the
        # protein of interest in the Protein Data Bank
        if do_get_structures:

            # Launch the process
            output_files = \
                client.submit(\
                    run_pipeline_step,
                    step = "get_structures",
                    wd = curr_wd,
                    uniprot_id = uniprot_id,
                    input_files = None,
                    output_files = [protein_file],
                    configrun_file = configrun_file,
                    log_file = LOG_GET_STRUCTURES_FILE,
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

            # If slimfast_find_slims will be run
            if do_find_slims:

                # Add the future to the dictionary storing the
                # futures slimfast_find_slims depends on
                tmp_futures["gs"] = output_files

            # Otherwise
            else:
                
                # Add the future to the list of futures that will
                # be gathered in the end
                futures.append(output_files)



        # If the user requested the search for structures of the
        # protein of interest in the AlphaFold Protein Structure
        # Database
        if do_get_alphafold:

            # Launch the process
            output_files = \
                client.submit(\
                    run_pipeline_step,
                    step = "get_alphafold",
                    wd = curr_wd,
                    uniprot_id = uniprot_id,
                    input_files = None,
                    output_files = [alphafold_file],
                    configrun_file = configrun_file,
                    log_file = LOG_GET_ALPHAFOLD_FILE,
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


        # If the user requested the search for structures of the
        # homologs of the protein of interest in the Protein Data
        # Bank
        if do_get_homologs:

            # Launch the process
            output_files = \
                client.submit(\
                    run_pipeline_step,
                    step = "get_homologs",
                    wd = curr_wd,
                    uniprot_id = uniprot_id,
                    input_files = None,
                    output_files = [homologs_file, alignments_file],
                    configrun_file = configrun_file,
                    log_file = LOG_GET_HOMOLOGS_FILE,
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


        # If the user requested the SLiMs' search
        if do_find_slims:

            # Get the input files
            input_files = \
                client.submit(\
                     util.get_items,
                        d = tmp_futures,
                        keys = ["gs", "ga", "gh"],
                        expected_num_values = [1, 1, 2],
                        insert_first = [slims_input_file])

            # Get the output files for slimfast_find_slims
            output_files = \
                [slims_protein_file, slims_alphafold_file,
                 slims_homologs_file, slims_disprot_file]

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
                    log_file = LOG_FIND_SLIMS_FILE,
                    v = v,
                    vv = vv)

            # Add the future to the list of futures that will
            # be gathered in the end
            futures.append(output_files)


    # Gather the remaining futures
    client.gather(futures)



def main():

    # Try to run
    try:

        run()

    # If something went wrong, report and exit
    except Exception as e:

        errstr = \
            f"Could not run slimfast_run_pipeline. Exception: {e}"
        logger.exception(errstr)
        sys.exit(errstr)

