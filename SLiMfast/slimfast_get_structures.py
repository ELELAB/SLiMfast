#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_get_structures.py
#
#    Given the UniProt ID of a protein, get data about its available
#    PDB structures and save such data in a CSV file.
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
import argparse
import logging as log
import os
import sys
# Third-party packages
import yaml
# SLiMfast
from SLiMfast import util



# Get the module logger
logger = log.getLogger(__name__)



def get_args():


    #------------------------ Argument parser ------------------------#

    
    # Set the script description
    description = \
        "\nGiven the UniProt ID of a protein, get data about its " \
        "available PDB structures and save such data in a CSV file.\n"

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
                              log_file_default = log_file_default)

    # Retrieve the arguments
    uniprot_id = args.uniprot_id
    output_protein_file = args.output_protein_file
    configrun_file = args.configrun_file
    log_file = args.log_file
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
    
    # Basic logging configuration
    log.basicConfig(level = level, 
                    filename = log_file, 
                    filemode = "w")


    #------------------------- Configuration -------------------------#


    # Try to load the configuration file for running
    try:
        config_run = yaml.safe_load(open(configrun_file, "r"))

    # If something went wrong, report it and exit   
    except Exception as e:
        errstr = f"Could not parse the run configuration file " \
                 f"{configrun_file}: {e}"
        logger.exception(errstr)
        sys.exit(errstr)


    # Return UniProt ID, output file and configuration
    return uniprot_id, output_protein_file, config_run



def run(uniprot_id, 
        output_protein_file,
        config_run):


    # Get the working directory where the pipeline will be run
    wd = util.get_wd(config_run)

    # Get the configuration for the PDB structures
    struct_params = config_run["structures"]["protein"]


    #------------------ Get the structures available -----------------#

    
    # Get all PDB structures (mapped to the correct chains)
    # containing the protein whose sequence is represented by
    # the UniProt ID
    try:
        structures_df = util.get_structures_data_protein(\
                            uniprot_id = uniprot_id,
                            **struct_params["df"])
    
    # If something went wrong, report it and exit
    except Exception as e:
        errstr = f"Could not obtain the PDB structures for " \
                 f"UniProt ID {uniprot_id}: {e}"
        logger.exception(errstr)
        sys.exit(errstr)

    # Specify the path to the output CSV file
    output_protein_path = os.path.join(wd, output_protein_file)

    # Save the dataframe to the output CSV file
    structures_df.to_csv(output_protein_path, **struct_params["out"])

    # Return the dataframe and the output CSV file
    return structures_df, output_protein_path



def main():
    args = get_args()
    run(*args)