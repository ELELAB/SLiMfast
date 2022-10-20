#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_get_structures.py
#
#    Given the UniProt ID of a protein, get data about its available
#    structures in the Protein Data Bank and save such data
#    in a CSV file.
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
from SLiMfast import util



# Get the module logger
logger = log.getLogger(__name__)



def get_args():


    #------------------------ Argument parser ------------------------#

    
    # Set the script description
    description = \
        "\nGiven the UniProt ID of a protein, get data about its " \
        "available structures in the Protein Data Bank and save " \
        "the data in a CSV file.\n"

    # Set the default name for the log file
    log_file_default = "slimfast_get_structures.log"

    # Add relevant arguments to the argument parser and return them
    args = util.get_arguments(description = description,
                              uniprot_id = True,
                              output_protein_file = True,
                              configrun_file = True,
                              log_file = True,
                              v = True,
                              vv = True,
                              log_to_stdout = True,
                              log_file_default = log_file_default)

    # Retrieve the arguments
    uniprot_id = args.uniprot_id
    output_protein_file = args.output_protein_file
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
    return uniprot_id, output_protein_file, config_run



def run(uniprot_id, 
        output_protein_file,
        config_run):


    # Get the working directory where the pipeline will be run
    wd = util.get_wd(config_run)

    # Get the configuration for the structures of the protein of
    # interest to be found in the Protein Data Bank
    struct_params = config_run["structures"]["protein"]


    #------------------ Get the structures available -----------------#

    
    # Get all structures (mapped to the correct chains) containing
    # the protein whose sequence is represented by the UniProt ID
    # in the Protein Data Bank
    try:
        
        structures_df = \
            util.get_structures_data_protein(\
                uniprot_id = uniprot_id,
                **struct_params["df"])
    
    # If something went wrong
    except Exception as e:

        # Report it and exit
        errstr = \
            f"Could not retrieve the data about the structures " \
            f"available in the Protein Data Bank for the protein " \
            f"with UniProt ID {uniprot_id}: {e}"
        logger.exception(errstr)
        sys.exit(errstr)

    # If no structures were found (fitting the criteria provided in
    # the configuration file)
    if structures_df.empty:
        
        # Inform the user and do not generate any output file
        warnstr = \
            f"No structures were found in the Protein Data Bank " \
            f"for the protein with UniProt ID {uniprot_id}."
        logger.warning(warnstr)

    # If structures were found
    else:

        # Specify the path to the output CSV file
        output_protein_path = os.path.join(wd, output_protein_file)

        # Save the data frame to the output CSV file
        structures_df.to_csv(output_protein_path,
                             **struct_params["out"])

        # Inform the user that the file was successfully written
        infostr = \
            f"The data frame containing data about the structures " \
            f"of the protein with UniProt ID {uniprot_id} found " \
            f"in the Protein Data Bank has been successfully " \
            f"written to {output_protein_path}."
        logger.info(infostr)

        # Return the data frame and the output CSV file
        return structures_df, output_protein_path



def main():
    args = get_args()
    run(*args)