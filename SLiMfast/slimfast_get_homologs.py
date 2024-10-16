#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_get_homologs.py
#
#    Given the UniProt ID of a protein, find its homologs, get data
#    about the available PDB structures of the homologs and save 
#    such data in a CSV file, and the alignments used to find the
#    homologs in a JSON file.
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
# SLiMfast
from .argparse_tools import get_arguments
from . import (
    defaults,
    dataframe_tools,
    hhsearch_tools,
    run_tools,
    util,
    )


#######################################################################


def run():

    # Set the script's name
    name = "slimfast_get_homologs"

    # Set the script's description
    description = \
        "\nGiven the UniProt ID of a protein, find its homologs, " \
        "get data about the structures of the homologs in the " \
        "in the Protein Data Bank and save the data in a CSV file. " \
        "Save the alignments used to find the homologs in a JSON " \
        "file.\n"

    # Set the names of the relevant arguments to be parsed from the
    # command line
    args_names = \
        ["uniprot_id", "output_homologs_file",
         "output_alignments_file", "configrun_file",
         "log_file", "v", "vv", "log_to_stdout"]

    # Set the default values for the arguments
    default_values = \
        {"log_file" : defaults.LOG_FILES["slimfast_get_homologs"]}

    #-----------------------------------------------------------------#

    # Add relevant arguments to the argument parser and return them
    args = get_arguments(name = name,
                         description = description,
                         args_names = args_names,
                         default_values = default_values)

    #-----------------------------------------------------------------#

    # Retrieve the arguments
    uniprot_id = args.uniprot_id
    output_homologs_file = args.output_homologs_file
    output_alignments_file = args.output_alignments_file
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
                             config_type = "slimfast_get_homologs")

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

    # Get the working directory where the program will be run
    wd = config_run["wd"]

    #-----------------------------------------------------------------#

    # Inform the user that the HHpred pipeline is about to run
    infostr = \
        "Preparing to run the HHpred pipeline (HHblits + HHsearch)..."
    log.info(infostr)

    #-----------------------------------------------------------------#

    # Try to run the HHpred pipeline
    try:
        
        _, hhsearch_output = \
            run_tools.run_hhpred(uniprot_id = uniprot_id,
                                 hhblits_params = \
                                    config_run["hhblits"]["run"],
                                 hhsearch_params = 
                                    config_run["hhsearch"]["run"],
                                 wd = wd)

    # If something went wrong
    except Exception as e:

        # Set the error message
        errstr = "It was not possible to complete the HHpred run."
        
        # Re-raise the exception
        raise Exception(errstr) from e

    # Inform the user that the HHpred run completed successfully
    infostr = "The HHpred run completed successfully."
    log.info(infostr)

    #-----------------------------------------------------------------#
    
    # Create the parser to parse the HHsearch output
    parser = hhsearch_tools.HHsearchParser()

    # Try to get the data frame with data about the available
    # structures and the alignments
    try:
        
        df, als = parser.parse(hhsearch_output = hhsearch_output)

    # If something went wrong
    except Exception as e:

        # Set the error message
        errstr = \
            "It was not possible to create the data frame " \
            "containing data about the structures of the homologs " \
            "of the protein of interest available in the Protein " \
            "Data Bank."
        
        # Re-raise the exception
        raise Exception(errstr) from e

    # Inform the user that the data frame was successfully created
    infostr = \
        "The data frame containing data about the structures of the " \
        "homologs of the protein of interest available in the " \
        "Protein Data Bank was successfully created."
    log.info(infostr)

    #-----------------------------------------------------------------#

    # Convert the index of the data frame into one of the columns
    # storing the unique IDs of the hits
    df = df.reset_index()
    df = df.rename(columns = {"index" : defaults.DF_COLUMNS["hit_id"]})

    #-----------------------------------------------------------------#

    # Try to add metadata to the data frame
    try:
        
        df = util.add_structures_metadata(df = df,
                                          source = "pdb")
    
    # If something went wrong
    except Exception as e:

        # Set the error message
        errstr = \
            "It was not possible to add metadata to the data frame."
        
        # Re-raise the exception
        raise Exception(errstr) from e

    # Inform the user that the metadata were successfully added
    infostr = "Metadata were successfully added to the data frame."
    log.info(infostr)

    #-----------------------------------------------------------------#

    # If there is a configuration for filtering the data frame
    if "filter" in config_run["structures"]["homolog"]["df"]:

        # Get the configuration
        config_filter = \
            config_run["structures"]["homolog"]["df"]["filter"]

        # Try to filter the data frame
        try:
            
            df = dataframe_tools.filter_df(df = df,
                                           filters = config_filter)

            # Only keep the hits in the filtered data frame in the
            # alignments
            als = \
                {hit : data for hit, data in als.items() \
                 if hit in df[defaults.DF_COLUMNS["hit_id"]].tolist()}

        # If something went wrong
        except Exception as e:

            # Set the error message
            errstr = "It was not possible to filter the data frame."

            # Re-raise the exception
            raise Exception(errstr) from e

        # Inform the user that the data frame was successfully
        # filtered
        infostr = "The data frame was successfully filtered."
        log.info(infostr)

    #-----------------------------------------------------------------#

    # If no structures were left after filtering
    if df.empty:

        # Inform the user and do not generate any output file
        warnstr = \
            "No structures were left in the data frame after " \
            "filtering."
        log.warning(warnstr)

        # Return
        return

    #-----------------------------------------------------------------#

    # If there is a configuration for sorting the data frame
    if "sort" in config_run["structures"]["homolog"]["df"]:

        # Get the configuration
        config_sort = \
            config_run["structures"]["homolog"]["df"]["sort"]

        # Try to sort the dataframe according to the given criteria
        try:
            
            df = \
                df.sort_values(by = config_sort["sort_by"], 
                               ascending = config_sort["ascending"])
        
        # If something went wrong
        except Exception as e:

            # Set the error message
            errstr = "It was not possible to sort the data frame."

            # Re-raise the exception
            raise Exception(errstr) from e

        # Inform the user that the data frame was successfully
        # sorted
        infostr = "The data frame was successfully sorted."
        log.info(infostr)

    #-----------------------------------------------------------------#

    # Specify the path to the output CSV file for the homologs'
    # data frame
    output_homologs_path = os.path.join(wd, output_homologs_file)

    # Save the data frame to the output CSV file
    df.to_csv(output_homologs_path,
              **defaults.OUTPUT_FILES["homologs"]["options"])

    # Inform the user that the file was successfully written
    infostr = \
        "The data frame was successfully saved in " \
        f"'{output_homologs_path}'."
    log.info(infostr)

    #-----------------------------------------------------------------#

    # Specify the path to the output JSON file for the alignments
    output_alignments_path = os.path.join(wd, output_alignments_file)

    # Save the alignments to the output JSON file
    with open(output_alignments_path, "w") as out:
        json.dump(als, out)

    # Inform the user that the file was successfully written
    infostr = \
        "The alignments of the structures of the homologs of " \
        "the protein of interest (aligned with the protein " \
        "itself) were successfully saved in " \
        f"'{output_alignments_path}'."
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
            f"'slimfast_get_homologs'. Error: {e}"
        log.exception(errstr)
        sys.exit(errstr)
