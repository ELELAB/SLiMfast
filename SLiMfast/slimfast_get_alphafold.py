#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_get_alphafold.py
#
#    Given the UniProt ID of a protein, get data about its available
#    structures in the AlphaFold Protein Structure Database
#     and save such data in a CSV file.
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
# SLiMfast
from . import dataframe_tools
from . import util
from .defaults import (
    ALPHAFOLD_ID_FIELD,
    DEPOSITION_DATE_FIELD,
    LOG_GET_ALPHAFOLD_FILE,
    STRUCT_CHAIN_FIELD,
    STRUCT_FRAGMENT_FIELD,
    STRUCT_ID_FIELD,
    STRUCT_SEGMENTS_FIELD,
    STRUCT_VERSION_FIELD,
    UNP_ID_FIELD,
    UNP_ISOFORM_FIELD,
    UNP_SEGMENTS_FIELD,
    )



# Get the module logger
logger = log.getLogger(__name__)



def get_args():


    #------------------------ Argument parser ------------------------#

    
    # Set the script description
    description = \
        "\nGiven the UniProt ID of a protein, get data about its " \
        "available structures in the AlphaFold Protein Structure " \
        "database and save the data in a CSV file.\n"

    # Add relevant arguments to the argument parser and return them
    args = \
        util.get_arguments(description = description,
                           uniprot_id = True,
                           output_alphafold_file = True,
                           configrun_file = True,
                           log_file = True,
                           v = True,
                           vv = True,
                           log_to_stdout = True,
                           log_file_default = LOG_GET_ALPHAFOLD_FILE)

    # Retrieve the arguments
    uniprot_id = args.uniprot_id
    output_alphafold_file = args.output_alphafold_file
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

    # If something went wrong 
    except Exception as e:
        
        # Report it and exit
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

    # Return UniProt ID, output file and configuration
    return uniprot_id, output_alphafold_file, config_run


def run(uniprot_id, 
        output_alphafold_file,
        config_run):


    # Get the working directory where the pipeline will be run
    wd = util.get_wd(config_run)

    # Get the configuration for the structures of the protein of
    # interest to be found in the AlphaFold Protein Structure Database
    struct_params = config_run["structures"]["alphafold"]


    #------------------ Get the structures available -----------------#


    # Sorted column names for the final data frame
    COLNAMES = \
        [STRUCT_ID_FIELD, STRUCT_CHAIN_FIELD, STRUCT_SEGMENTS_FIELD,
         STRUCT_VERSION_FIELD, STRUCT_FRAGMENT_FIELD,
         UNP_ID_FIELD, UNP_SEGMENTS_FIELD, UNP_ISOFORM_FIELD, 
         ALPHAFOLD_ID_FIELD, DEPOSITION_DATE_FIELD]


    #------------------ Get the available structures -----------------#


    # Try to get the data frame of available structures
    try:
        
        df = util.get_alphafold_df(uniprot_id = uniprot_id)

        # Inform the user that the IDs of the available structures
        # of the protein were successfully retrieved
        infostr = \
            f"The IDs of the available structures of the protein " \
            f"with UniProt ID {uniprot_id} have been successfully " \
            f"retrieved from the AlphaFold Protein Structure Database."
        logger.info(infostr)
    
    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not get the IDs of the structures available in " \
            f"the AlphaFold Protein Structure Database for the " \
            f"protein with UniProt ID {uniprot_id}. Exception: {e}"
        logger.exception(errstr)

        # Re-raise the exception
        raise Exception(errstr)

    # If no structures were found
    if df.empty:

        # Inform the user and do not generate any output file
        warnstr = \
            f"No structures were found in the AlphaFold Protein " \
            f"Structure Database for the protein with UniProt ID " \
            f"{uniprot_id}."
        logger.warning(warnstr)

        # Return
        return None, None


    #------------------------ Add the metadata -----------------------#


    # Try to add metadata for the structures
    try:
        
        df = util.add_structures_metadata(df = df,
                                          source = "alphafold")
    
    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not add metadata to the structures found for " \
            f"the protein with UniProt ID {uniprot_id} in the " \
            f"AlphaFold Protein Structure Database. Exception: {e}"
        logger.exception(errstr)

        # Re-raise the exception
        raise Exception(errstr)


    #---------------------- Filter the data frame --------------------#


    # Try to filter the data frame
    try:
        
        df = dataframe_tools.filter_df(\
                df = df,
                filters = struct_params["df"]["filter"])

        # Inform the user that the data frame was successfully
        # filtered
        infostr = \
            f"The data frame of the structures of the protein " \
            f"with UniProt ID {uniprot_id} found in the AlphaFold " \
            f"Protein Structure Database has been successfully " \
            f"filtered."
        logger.info(infostr)

    # If something went wrong
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not filter the data frame of the structures " \
            f"of the protein with UniProt ID {uniprot_id} found in " \
            f"the AlphaFold Protein Structure Database. Exception: {e}"
        logger.exception(errstr)

        # Re-raise the exception
        raise Exception(errstr)


    # If no structures remained after filtering
    if df.empty:

        # Inform the user and do not generate any output file
        warnstr = \
            f"No structures were found in the AlphaFold Protein " \
            f"Structure Database for the protein with UniProt ID " \
            f"{uniprot_id}."
        logger.warning(warnstr)

        # Return
        return None, None


    #---------------------- Sort the data frame ----------------------#


    # Try to sort the dataframe according to the given criteria
    try:

        df = df.sort_values(\
            by = struct_params["df"]["sort"]["sort_by"], 
            ascending = struct_params["df"]["sort"]["ascending"])

        # Inform the user that the data frame was successfully
        # sorted
        infostr = \
            f"The data frame of the structures of the protein " \
            f"with UniProt ID {uniprot_id} found in the AlphaFold " \
            f"Protein Structure Database has been successfully sorted."
        logger.info(infostr)
 
     # If something went wrong 
    except Exception as e:

        # Log the error
        errstr = \
            f"Could not sort the data frame of the structures of " \
            f"the protein with UniProt ID {uniprot_id} found in " \
            f"the AlphaFold Protein Structure Database. " \
            f"Exception: {e}"
        logger.exception(errstr)

        # Re-raise the exception
        raise Exception(errstr)

    # Inform the user that all data were successfully collected
    infostr = \
        f"All data for the structures of the protein with " \
        f"UniProt ID {uniprot_id} found in the AlphaFold Protein " \
        f"Structure Database have been successfully collected."
    logger.info(infostr)


    #----------------------- Sort the columns ------------------------#


    # Sort the columns
    df = df.reindex(columns = COLNAMES)


    #----------------------- Write the CSV file ----------------------#


    # Specify the path to the output CSV file
    output_alphafold_path = os.path.join(wd, output_alphafold_file)

    # Save the data frame to the output CSV file
    df.to_csv(output_alphafold_path,
              **struct_params["out"])

    # Inform the user that the file was successfully written
    infostr = \
        f"The data frame containing data about the structures " \
        f"of the protein with UniProt ID {uniprot_id} found " \
        f"in the AlphaFold Protein Structure Database has been " \
        f"successfully written to {output_alphafold_path}."
    logger.info(infostr)

    # Return the data frame and the output CSV file
    return df, output_alphafold_path



def main():
    
    # Try to run
    try:

        args = get_args()
        run(*args)

    # If something went wrong, report and exit
    except Exception as e:

        errstr = \
            f"Could not run slimfast_get_alphafold. Exception: {e}"
        logger.exception(errstr)
        sys.exit(errstr)
