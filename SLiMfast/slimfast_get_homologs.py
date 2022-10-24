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
from SLiMfast import util



# Get the module logger
logger = log.getLogger(__name__)



def get_args():


    #------------------------ Argument parser ------------------------#
    

    # Set the script description
    description = \
        "\nGiven the UniProt ID of a protein, find its homologs, " \
        "get data about the structures of the homologs in the " \
        "in the Protein Data Bank and save the data in a CSV file. " \
        "Save the alignments used to find the homologs in a JSON " \
        "file.\n"

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
    # resulting data frame and alignments
    try:
        structures_df, als = \
            util.get_structures_data_homologs(\
                uniprot_id = uniprot_id,
                hhblits_params = hhblits_params,
                hhsearch_params = hhsearch_params,
                wd = wd,
                **struct_params["df"])
    
    # If something went wrong
    except Exception as e:
        
        # Report it and exit
        errstr = \
            f"Could not retrieve the data about the structures " \
            f"available in the Protein Data Bank for the homologs " \
            f"of the protein with UniProt ID {uniprot_id}: {e}"
        logger.exception(errstr)
        sys.exit(errstr)

    # If no homologs were found (fitting the criteria provided in
    # the configuration file)
    if structures_df.empty:

        # Inform the user and do not generate any output file
        warnstr = \
            f"No structures were found in the Protein Data Bank " \
            f"matching the criteria provided in the configuration " \
            f"file for the protein with UniProt ID {uniprot_id}. " \
            f"Relaxing the filtering criteria may help in getting " \
            f"more hits."
        logger.warning(warnstr)

    # If homologs were found
    else:

        # Specify the path to the output CSV file for the homologs'
        # data frame and save it to the CSV file
        output_homologs_path = \
            os.path.join(wd, output_homologs_file)
        structures_df.to_csv(output_homologs_path,
                             **struct_params["out"])

        # Inform the user that the file was successfully written
        infostr = \
            f"The data frame containing data about the structures " \
            f"of the homologs of the protein with UniProt ID " \
            f"{uniprot_id} found in the Protein Data Bank has been " \
            f"successfully written to {output_homologs_path}."
        logger.info(infostr)

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

        # Return the homologs' data frame, the alignments and the paths
        # to the output files
        return (structures_df, als), \
               (output_homologs_path, output_alignments_path)



def main():
    args = get_args()
    run(*args)