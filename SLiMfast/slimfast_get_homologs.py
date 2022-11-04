#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_get_homologs.py
#
#    Given the UniProt ID of a protein, find its homologs, get data
#    about the available PDB structures of the homologs and save 
#    such data in a CSV file, and the alignments used to find the
#    homologs in a JSON file.
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
# SLiMfast
from .argparse_tools import get_arguments
from . import (
    dataframe_tools,
    hhsearch_tools,
    run_tools,
    util,
    )
from .defaults import (
    HIT_ID_FIELD,
    LOG_GET_HOMOLOGS_FILE,
    )



# Get the module logger
logger = log.getLogger(__name__)



def run():



    #------------------------ Argument parser ------------------------#



    # Set the script description
    description = \
        "\nGiven the UniProt ID of a protein, find its homologs, " \
        "get data about the structures of the homologs in the " \
        "in the Protein Data Bank and save the data in a CSV file. " \
        "Save the alignments used to find the homologs in a JSON " \
        "file.\n"

    # Add relevant arguments to the argument parser and return them
    args = \
        get_arguments(\
            description = description,
            uniprot_id = True,
            output_homologs_file = True,
            output_alignments_file = True,
            configrun_file = True,
            log_file = True,
            v = True,
            vv = True,
            log_to_stdout = True,
            log_file_default = LOG_GET_HOMOLOGS_FILE)

    # Retrieve the arguments
    uniprot_id = args.uniprot_id
    output_homologs_file = args.output_homologs_file
    output_alignments_file = args.output_alignments_file
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
    
    # If something went wrong, report it and exit
    except Exception as e:
        errstr = \
            f"Could not parse the configuration file " \
            f"{configrun_file}: {e}"
        logger.exception(errstr)
        sys.exit(errstr)

    # Inform the user about the successful parsing of the
    # configuration file
    infostr = \
        f"Configuration parsed successfully from {configrun_file}."
    logger.info(infostr)



    #------------------------------ Run ------------------------------#



    # Get the working directory where the pipeline will be run
    wd = util.get_wd(config_run)

    # Get the configuration for the PDB structures
    struct_params = config_run["structures"]["homolog"]

    # Get the HHblits and HHsearch parameters
    hhblits_params = config_run["hhblits"]["run"]
    hhsearch_params = config_run["hhsearch"]["run"]



    #--------------------------- Run HHpred --------------------------#



    # Inform the user that the HHpred pipeline is about to run
    infostr = \
        "Preparing to run the HHpred pipeline (HHpred + HHsearch)..."
    logger.info(infostr)

    # Try to run the HHpred pipeline
    try:
        
        _, hhsearch_output = \
            run_tools.run_hhpred(uniprot_id = uniprot_id,
                                 hhblits_params = hhblits_params,
                                 hhsearch_params = hhsearch_params,
                                 wd = wd)

    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = f"Could not complete the HHpred run. Exception: {e}"
        logger.exception(errstr)

        # Re-raise the exception
        raise Exception(errstr)



    #------------------ Get the homologs' structures -----------------#



    # Create the parser to parse the HHsearch output
    parser = hhsearch_tools.HHsearchParser()

    # Try to parse the HHsearch output to get the list of homologs
    # of the protein of interest and their associated data
    try:
        
        df, als = parser.parse(hhsearch_output = hhsearch_output)

        # Inform the user that the IDs of the available structures
        # for the homologs were successfully retrieved
        infostr = \
            f"The IDs of the structures of the homologs of the " \
            f"protein with UniProt ID {uniprot_id} in the Protein " \
            f"Data Bank have been successfully retrieved."
        logger.info(infostr)

    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not get the structures of the homologs of the " \
            f"protein with UniProt ID {uniprot_id} found in the " \
            f"Protein Data Bank. Exception: {e}"
        logger.exception(errstr)
        
        # Re-raise the exception
        raise Exception(errstr)

    # Convert the index of the data frame into one of the columns
    # storing the unique IDs of the hits
    df = df.reset_index()
    df = df.rename(columns = {"index" : HIT_ID_FIELD})



    #------------------------ Add the metadata -----------------------#



    # Try to add structures' metadata
    try:
        
        df = util.add_structures_metadata(df = df,
                                          source = "pdb")
    
    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not add metadata to the structures of the " \
            f"homologs of the protein with UniProt ID {uniprot_id} " \
            f"found in the Protein Data Bank. Exception: {e}"
        logger.exception(errstr)
        
        # Re-raise the exception
        raise Exception(errstr)



    #---------------------- Filter the data frame --------------------#



    # Try to filter the hits
    try:
        
        df = dataframe_tools.filter_df(\
                df = df,
                filters = struct_params["df"]["filter"])

        # Only keep the hits in the filtered data frame in the
        # alignments
        als = \
            {hit : data for hit, data in als.items() \
             if hit in df[HIT_ID_FIELD].tolist()}

        # Inform the user that the data frame was successfully
        # filtered
        infostr = \
            f"The data frame of the structures of the homologs of " \
            f"protein with UniProt ID {uniprot_id} found in the " \
            f"Protein Data Bank has been successfully filtered."
        logger.info(infostr)
    
    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not filter the data frame of the structures " \
            f"of the homologs of the protein with UniProt ID " \
            f"{uniprot_id} found in the Protein Data Bank. " \
            f"Exception: {e}"
        logger.exception(errstr)
        
        # Re-raise the exception
        raise Exception(errstr)



    #----------------------- Sort the datarame -----------------------#



    # Try to sort the hits
    try:
        
        df = \
            df.sort_values(\
                by = struct_params["df"]["sort"]["sort_by"], 
                ascending = struct_params["df"]["sort"]["ascending"])

        # Inform the user that the data frame was successfully
        # sorted
        infostr = \
            f"The data frame of the structures of the homologs of " \
            f"the protein with UniProt ID {uniprot_id} found in the " \
            f"Protein Data Bank has been successfully sorted."
        logger.info(infostr)
    
    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not sort the data frame of the structures of " \
            f"the homologs of the protein with UniProt ID " \
            f"{uniprot_id} found in the Protein Data Bank. " \
            f"Exception: {e}"
        logger.exception(errstr)
        
        # Re-raise the exception
        raise Exception(errstr)

    infostr = \
        f"All data for the structures the homologs of the protein " \
        f"with {uniprot_id} found in the Protein Data Bank " \
        f"have been successfully collected."
    logger.info(infostr)


    # If no homologs were found (fitting the criteria provided in
    # the configuration file)
    if df.empty:

        # Inform the user and do not generate any output file
        warnstr = \
            f"{uniprot_id}: no structures remained in the data " \
            f"frame after filtering. Relaxing the filtering " \
            f"criteria may help in getting more hits."
        logger.warning(warnstr)

        # Return
        return



    #----------------------- Write the CSV file ----------------------#



    # Specify the path to the output CSV file for the homologs'
    # data frame and save it to the CSV file
    output_homologs_path = \
        os.path.join(wd, output_homologs_file)
    df.to_csv(output_homologs_path,
              **struct_params["out"])

    # Inform the user that the file was successfully written
    infostr = \
        f"The data frame containing data about the structures " \
        f"of the homologs of the protein with UniProt ID " \
        f"{uniprot_id} found in the Protein Data Bank has been " \
        f"successfully written to {output_homologs_path}."
    logger.info(infostr)



    #---------------------- Write the JSON file ----------------------#



    # Specify the path to the output JSON file for the alignments
    # and save them to the JSON file
    output_alignments_path = \
        os.path.join(wd, output_alignments_file)
    with open(output_alignments_path, "w") as out:
        json.dump(als, out)

    # Inform the user that the file was successfully written
    infostr = \
        f"The alignments of the structures of the homologs of " \
        f"the protein with UniProt ID {uniprot_id} " \
        f"(aligned with the protein itself) found in the " \
        f"Protein Data Bank have been successfully written " \
        f"to {output_alignments_path}."
    logger.info(infostr)



def main():

    # Try to run
    try:

        run()

    # If something went wrong, report and exit
    except Exception as e:

        errstr = \
            f"Could not run slimfast_get_homologs. Exception: {e}"
        logger.exception(errstr)
        sys.exit(errstr)