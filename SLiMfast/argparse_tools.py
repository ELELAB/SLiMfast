#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    argparse_tools.py
#
#    Utilities to handle the command-line arguments passed to the
#    SLiMfast executables.
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


__doc__ = \
    "Utilities to handle the command-line arguments passed to " \
    "the SLiMfast executables."


#######################################################################


# Standard library
import argparse
# SLiMfast
from . import defaults


#######################################################################


def get_arguments(name,
                  description,
                  args_names,
                  default_values = None):
    """Create an argument parser for a script, add the requested
    arguments, and return them.

    Parameters
    ----------
    name : ``str``
        The name of the script.

    description : ``str``
        The description of the script.

    args_names : ``list``
        A list containing the names of the arguments to be added
        to the parser.

    default_values : ``dict``, optional
        A dictionary mapping arguments' names to their default
        values, if any.

    Returns
    -------
    args : ``argparser.Namespace``
        The namespace where the arguments live.
    """

    # Create the argument parser
    parser = argparse.ArgumentParser(prog = name,
                                     description = description)

    #-----------------------------------------------------------------#
    
    # The UniProt ID of the protein of interest.
    #
    # Taken by:
    # - slimfast_get_structures
    # - slimfast_get_homologs
    if "uniprot_id" in args_names:

        # Set the help message
        help_msg = "The UniProt ID of the protein of interest."

        # Add the argument to the parser
        parser.add_argument("-id", "--uniprot-id",
                             type = str,
                             required = True,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # Whether to run 'slimfast_get_structures' as part of the pipeline.
    #
    # Taken by:
    # - slimfast_run_pipeline
    if "do_get_structures" in args_names:

        # Set the help message
        help_msg = "Run 'slimfast_get_structures'."

        # Add the argument to the parser
        parser.add_argument("-dgs", "--do-get-structures",
                            action = "store_true",
                            help = help_msg)

    #-----------------------------------------------------------------#

    # Whether to run 'slimfast_get_alphafold' as part of the pipeline.
    #
    # Taken by:
    # - slimfast_run_pipeline
    if "do_get_alphafold" in args_names:

        # Set the help message
        help_msg = "Run 'slimfast_get_alphafold'."

        # Add the argument to the parser
        parser.add_argument("-dga", "--do-get-alphafold",
                            action = "store_true",
                            help = help_msg)

    #-----------------------------------------------------------------#

    # Whether to run 'slimfast_get_homologs' as part of the pipeline.
    #
    # Taken by:
    # - slimfast_run_pipeline
    if "do_get_homologs" in args_names:

        # Set the help message
        help_msg = "Run 'slimfast_get_homologs'."

        # Add the argument to the parser
        parser.add_argument("-dgh", "--do-get-homologs",
                            action = "store_true",
                            help = help_msg)

    #-----------------------------------------------------------------#

    # Whether to run 'slimfast_find_slims' as part of the pipeline.
    # 
    # Taken by:
    # - slimfast_run_pipeline
    if "do_find_slims" in args_names:

        # Set the help message
        help_msg = "Run 'slimfast_find_slims'."

        # Add the argument to the parser
        parser.add_argument("-dfs", "--do-find-slims",
                            action = "store_true",
                            help = help_msg)

    #-----------------------------------------------------------------#

    # The CSV file containing the UniProt IDs of the proteins of
    # interest and the associated SLiMs to be found in their
    # sequences.
    #
    # Taken by:
    # - slimfast_run_pipeline
    if "input_ids_slims_file" in args_names:

        # Set the help message
        help_msg = \
            "The input CSV file containing the UniProt IDs of the " \
            "proteins of interest and the SLiMs to be found " \
            "in their sequences."

        # Add the argument to the parser
        parser.add_argument("-ids", "--input-ids-slims-file",
                             type = str,
                             required = True,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The CSV file containing the structures of the protein of interest
    # found in the Protein Data Bank.
    #
    # Taken by:
    # - slimfast_find_slims
    #
    # Generated by:
    # - slimfast_get_structures
    if "input_protein_file" in args_names:

        # Set the help message
        help_msg = \
            "The input CSV file containing the structures of the " \
            "protein of interest found in the Protein Data Bank."

        # Add the argument to the parser
        parser.add_argument("-ip", "--input-protein-file",
                             type = str,
                             default = None,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The CSV file containing the structures of the protein of interest
    # found in the AlphaFold Protein Structure Database.
    #
    # Taken by:
    # - slimfast_find_slims
    #
    # Generated by:
    # - slimfast_get_alphafold
    if "input_alphafold_file" in args_names:

        # Set the help message
        help_msg = \
            "The input CSV file containing the structures of the " \
            "protein of interest found in the AlphaFold Protein " \
            "Structure Database."

        # Add the argument to the parser
        parser.add_argument("-ia", "--input-alphafold-file",
                             type = str,
                             default = None,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The CSV file containing the structures of the homologs of the
    # protein of interest found in the Protein Data Bank.
    #
    # Taken by:
    # - slimfast_find_slims
    #
    # Generated by:
    # - slimfast_get_homologs
    if "input_homologs_file" in args_names:

        # Set the help message
        help_msg = \
            "The input CSV file containing the structures of the " \
            "homologs of the protein of interest found in the " \
            "Protein Data Bank."

        # Add the argument to the parser
        parser.add_argument("-ih", "--input-homologs-file",
                             type = str,
                             default = None,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The JSON file containing the alignments between the homologs'
    # structures and the protein of interest.
    #
    # Taken by:
    # - slimfast_find_slims
    #    
    # Generated by:
    # - slimfast_get_homologs
    if "input_alignments_file" in args_names:

        # Set the help message
        help_msg = \
            "The input JSON file containing the alignments between " \
            "each of the structures of the homologs found in the " \
            "Protein Data Bank and the protein of interest."

        # Add the argument to the parser
        parser.add_argument("-ial", "--input-alignments-file",
                             type = str,
                             default = None,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The input CSV file containing the SLiMs found in the UniProt
    # sequence of the protien of interest.
    #
    # Taken by:
    # - slimfast_find_slims
    if "input_slims_file" in args_names:

        # Set the help message
        help_msg = \
            "The input CSV file containing the SLiMs in the UniProt " \
            "sequence of the protein of interest."

        # Add the argument to the parser
        parser.add_argument("-is", "--input-slims-file",
                            type = str,
                            required = True,
                            help = help_msg)

    #-----------------------------------------------------------------#

    # The name of the output CSV file containing the structures of the
    # protein of interest found in the Protein Data Bank.
    #
    # Taken by:
    # - slimfast_get_structures
    if "output_protein_file" in args_names:

        # Set the default value
        default = defaults.OUTPUT_FILES["protein"]["name"]

        # Set the help message
        help_msg = \
            "The name of the output CSV file containing the " \
            "structures of the protein of interest found in the " \
            "Protein Data Bank. The file is written in the 'wd' " \
            "directory specified in the configuration file. " \
            f"The default file name is '{default}'."

        # Add the argument to the parser
        parser.add_argument("-op", "--output-protein-file",
                             type = str,
                             default = default,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The name of the output CSV file containing the structures of the
    # protein of interest found in the AlphaFold Protein Structure
    # Database.
    #
    # Taken by:
    # - slimfast_get_alphafold
    if "output_alphafold_file" in args_names:

        # Set the default value
        default = defaults.OUTPUT_FILES["alphafold"]["name"]

        # Set the help message
        help_msg = \
            "The name of the output CSV file containing the " \
            "structures of the protein of interest found in the " \
            "AlphaFold Protein Structure Database. The file is " \
            "written in the 'wd' directory specified in the " \
            "configuration file. The default file name is " \
            f"'{default}'."

        # Add the argument to the parser
        parser.add_argument("-oa", "--output-alphafold-file",
                             type = str,
                             default = default,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The name of the output CSV file containing the structures of the
    # homologs of the protein of interest found in the Protein Data
    # Bank.
    #
    # Taken by:
    # - slimfast_get_homologs
    if "output_homologs_file" in args_names:

        # Set the default value
        default = defaults.OUTPUT_FILES["homologs"]["name"]

        # Set the help message
        help_msg = \
            "The name of the output CSV file containing the " \
            "structures of the homologs of the protein of " \
            "interest found in the Protein Data Bank. The file " \
            "is written in the 'wd' directory specified in the " \
            "configuration file. The default file name is " \
            f"'{default}'."

        # Add the argument to the parser
        parser.add_argument("-oh", "--output-homologs-file",
                             type = str,
                             default = default,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The name of the output JSON file containing the alignments
    # between the homologs' structures and the protein of interest.
    #
    # Taken by:
    # - slimfast_get_homologs
    if "output_alignments_file" in args_names:

        # Set the default value
        default = defaults.OUTPUT_FILES["alignments"]["name"]

        # Set the help message
        help_msg = \
            "The name of the output JSON file containing the " \
            "alignments between each of the structures of the " \
            "homologs found in the Protein Data Bank and the " \
            "protein of interest. The file is written in the 'wd' " \
            "directory specified in the configuration file. The " \
            f"default file name is '{default}'."

        # Add the argument to the parser
        parser.add_argument("-oal", "--output-alignments-file",
                             type = str,
                             default = default,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The name of the output CSV file containing the SLiMs found in
    # the structures of the protein of interest found in the Protein
    # Data Bank.
    #
    # Taken by:
    # - slimfast_find_slims
    if "output_slims_protein_file" in args_names:

        # Set the default value
        default = defaults.OUTPUT_FILES["slims_protein"]["name"]

        # Set the help message
        help_msg = \
            "The name of the output CSV file containing the SLiMs " \
            "found in the structures of the protein of interest " \
            "found in the Protein Data Bank. The file is written " \
            "in the 'wd' directory specified in the configuration " \
            f"file. The default file name is '{default}'."

        # Add the argument to the parser
        parser.add_argument("-osp", "--output-slims-protein-file",
                             type = str,
                             default = default,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The name of the output CSV file containing the SLiMs found in
    # the structures of the protein of interest found in the AlphaFold
    # Protein Structure Database.
    #
    # Taken by:
    # - slimfast_find_slims
    if "output_slims_alphafold_file" in args_names:

        # Set the default value
        default = defaults.OUTPUT_FILES["slims_alphafold"]["name"]

        # Set the help message
        help_msg = \
            "The name of the output CSV file containing the SLiMs " \
            "found in the structures of the protein of interest " \
            "found in the AlphaFold Protein Structure Database. " \
            "The file is written in the 'wd' directory specified " \
            "in the configuration file. The default file name is " \
            f"'{default}'."

        # Add the argument to the parser
        parser.add_argument("-osa", "--output-slims-alphafold-file",
                             type = str,
                             default = default,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The name of the output CSV file containing the SLiMs found in
    # the structures of the homolgs of the protein of interest found
    # in the Protein Data Bank.
    #
    # Taken by:
    # - slimfast_find_slims
    if "output_slims_homologs_file" in args_names:

        # Set the default value
        default = defaults.OUTPUT_FILES["slims_homologs"]["name"]

        # Set the help message
        help_msg = \
            "The name of the output CSV file containing the SLiMs " \
            "found in the structures of the homologs of the " \
            "protein of interest found in the Protein Data " \
            "Bank. The file is written in the 'wd' directory " \
            "specified in the configuration file. The default " \
            f"file name is '{default}'."

        # Add the argument to the parser
        parser.add_argument("-osh", "--output-slims-homologs-file",
                             type = str,
                             default = default,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The name of the output CSV file containing the SLiMs not found
    # in the structures and annotated through DisProt.
    #
    # Taken by:
    # - slimfast_find_slims
    if "output_slims_disprot_file" in args_names:

        # Set the default value
        default = defaults.OUTPUT_FILES["slims_disprot"]["name"]

        # Set the help message
        help_msg = \
            "The name of the output CSV file containing the SLiMs " \
            "not found in the structures available either for the " \
            "protein of interest or for the homologs, but annotated " \
            "as disordered in DisProt. The file is written in the " \
            "'wd' directory specified in the configuration file. " \
            f"The default file name is '{default}'."

        # Add the argument to the parser
        parser.add_argument("-osd", "--output-slims-disprot-file",
                             type = str,
                             default = default,
                             help = help_msg)

    #-----------------------------------------------------------------#

    # The YAML configuration file.
    #
    # Taken by:
    # - slimfast_get_structures
    # - slimfast_get_alphafold
    # - slimfast_get_homologs
    # - slimfast_find_slims
    if "configrun_file" in args_names:

        # Set the directory where the default files can be found
        default_dir = defaults.CONFIG_RUN_DIR

        # Set the help message
        help_msg = \
            "The YAML configuration file for running. If it is a " \
            "name without an extension, it is assumed to be a " \
            f"YAML file in '{default_dir}' with that name."

        # Add the argument to the parser
        parser.add_argument("-cr", "--configrun-file",
                            type = str,
                            required = True,
                            help = help_msg)

    #-----------------------------------------------------------------#

    # The name of the log file.
    #
    # Taken by:
    # - slimfast_get_structures
    # - slimfast_get_alphafold
    # - slimfast_get_homologs
    # - slimfast_find_slims
    if "log_file" in args_names:

        # Set the default value
        default = default_values.get("log_file")

        # Set the help message
        help_msg = \
            "The name of the log file. The file is written in the " \
            "'wd' directory specified in the configuration file. " \
            f"The default file name is '{default}'."

        # Add the argument to the parser
        parser.add_argument("-lf", "--log-file",
                            type = str,
                            default = default,
                            help = help_msg)

    #-----------------------------------------------------------------#

    # Whether to log also to the standard output.
    #
    # Taken by:
    # - slimfast_get_structures
    # - slimfast_get_alphafold
    # - slimfast_get_homologs
    # - slimfast_find_slims
    if "log_to_stdout" in args_names:

        # Set the help message
        help_msg = "Log also to the standard output."

        # Add the argument to the parser
        parser.add_argument("-ls", "--log-to-stdout",
                            action = "store_true",
                            help = help_msg)

    #-----------------------------------------------------------------#
    
    # Whether to enable verbose (INFO level) logging.
    #
    # Taken by:
    # - slimfast_get_structures
    # - slimfast_get_alphafold
    # - slimfast_get_homologs
    # - slimfast_find_slims
    if "v" in args_names:

        # Set the help message
        help_msg = "Whether to enable verbose (INFO level) logging."

        # Add the argument to the parser
        parser.add_argument("-v", "--verbose",
                            action = "store_true",
                            help = help_msg)

    #-----------------------------------------------------------------#

    # Whether to enable maximally verbose (DEBUG level) logging.
    #
    # Taken by:
    # - slimfast_get_structures
    # - slimfast_get_alphafold
    # - slimfast_get_homologs
    # - slimfast_find_slims
    if "vv" in args_names:

        # Set the help message
        help_msg = \
            "Whether to enable maximally verbose (DEBUG level) " \
            "logging."

        # Add the argument to the parser
        parser.add_argument("-vv", "--verbose-debug",
                            action = "store_true",
                            help = help_msg)

    #-----------------------------------------------------------------#

    # The number of processes to use when running.
    #
    # Taken by:
    # - slimfast_run_pipeline
    if "nproc" in args_names:

        # Set the default value
        default = defaults.NPROC_DEFAULT

        # Set the help message
        help_msg = \
            "The number of processes (cores) to use when running. " \
            f"The default is {default} process(es)."

        # Add the argument to the parser
        parser.add_argument("-n", "--nproc",
                            type = int,
                            default = default,
                            help = help_msg)

    #-----------------------------------------------------------------#

    # Parse the command line arguments
    args = parser.parse_args()

    #-----------------------------------------------------------------#

    # Return the arguments
    return args
