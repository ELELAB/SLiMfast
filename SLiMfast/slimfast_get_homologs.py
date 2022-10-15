#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    slimfast_get_homologs.py
#
#    Given the UniProt ID of a protein, find its homologs, get data
#    about the available PDB structures of the homologs and save 
#    such data in a CSV file, and the alignments used to find the
#    homologs in a JSON file.
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
import json
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
        "\nGiven the UniProt ID of a protein, find its homologs, " \
        "get data about the available PDB structures of the " \
        "homologs and save such data in a CSV file, and the " \
        "alignments used to find the homologs in a JSON file.\n"

    # Set the default name of the log file
    log_file_default = "slimfast_get_homologs.log"

    # Add relevant arguments to the argument parser and return them
    args = util.get_arguments(description = description,
                              uniprot_id = True,
                              output_homologs_file = True,
                              output_alignments_file = True,
                              configrun_file = True,
                              log_file = True,
                              v = True,
                              vv = True,
                              log_to_stdout = True,
                              log_file_default = log_file_default)

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
        config_run = yaml.safe_load(open(configrun_file, "r"))
    except Exception as e:
        # If something went wrong, report it and exit
        errstr = f"Could not parse the run configuration file " \
                 f"{configrun_file}: {e}"
        logger.exception(errstr)
        sys.exit(errstr)


    # Return UniProt ID, output files and configuration
    return uniprot_id, output_homologs_file, \
           output_alignments_file, config_run


def run(uniprot_id,
        output_homologs_file,
        output_alignments_file,
        config_run):


    # Get the working directory where the pipeline will be run
    wd = util.get_wd(config_run)

    # Get the configuration for the PDB structures
    struct_params = config_run["structures"]["homolog"]

    # Get the HHblits and HHsearch parameters
    hhblits_params = config_run["hhblits"]["run"]
    hhsearch_params = config_run["hhsearch"]["run"]


    #---------------------- Search for homologs ----------------------#


    # Run the HHpred pipeline (HHblits + HHsearch) and get the
    # resulting dataframe and alignments
    try:
        structures_df, als = \
            util.get_structures_data_homologs(\
                uniprot_id = uniprot_id,
                hhblits_params = hhblits_params,
                hhsearch_params = hhsearch_params,
                wd = wd,
                **struct_params["df"])
    
    # If something went wrong, report it and exit
    except Exception as e:
        errstr = f"Could not find homologs for the protein " \
                 f"of interest: {e}"
        logger.exception(errstr)
        sys.exit(errstr)

    # If no homologs were found (fitting the criteria provided in
    # the configuration file)
    if structures_df.empty:

        # Inform the user and do not generate any output file
        warnstr = \
            f"No structures were found matching the criteria " \
            f"provided in the configuration file. Relaxing the " \
            f"filtering criteria may help in getting more hits."
        logger.warning(warnstr)

    # If homologs were found
    else:

        # Specify the path to the output CSV file for the homologs'
        # dataframe and save it to the CSV file
        output_homologs_path = \
            os.path.join(wd, output_homologs_file)
        structures_df.to_csv(output_homologs_path,
                             **struct_params["out"])

        # Specify the path to the output JSON file for the alignments
        # and save them to the JSON file
        output_alignments_path = \
            os.path.join(wd, output_alignments_file)
        with open(output_alignments_path, "w") as out:
            json.dump(als, out)

        # Return the homologs' dataframe, the alignments and the paths
        # to the output files
        return (structures_df, als), \
               (output_homologs_path, output_alignments_path)



def main():
    args = get_args()
    run(*args)