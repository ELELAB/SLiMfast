#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    defaults.py
#
#    (Mostly internal) fefault values for SLiMfast.
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


# The module's description
__doc__ = "(Mostly internal) fefault values for SLiMfast."


#######################################################################


# Standard library
import logging as log
import os


#######################################################################


# Templates agsinst which to check the configurations
CONFIG_TEMPLATES = \
    {# The template for the configuration for 'slimfast_get_structures'
     "slimfast_get_structures" : \

        {# The directory where the executable will be run
         "wd" :  (str, os.getcwd()),
         
         # The settings to be used for the structures
         "structures" : \
            {# The settings to be used when running the executable
             "protein" : \
                {# The settings for the output data frame
                 "df" : \
                    {# The settings for sorting the data frame
                     "sort" : \
                        {# The sorting criteria
                         "sort_by" : list,
                         # For each sorting criterion, whether
                         # to sortin ascending order
                         "ascending" : list}}}},
        },

    #-----------------------------------------------------------------#

    # The template for the configuration for 'slimfast_get_alphafold'
     "slimfast_get_alphafold" : \

        {# The directory where the executable will be run
         "wd" :  (str, os.getcwd()),
         
         # The settings to be used for the structures
         "structures" : \
            {# The settings to be used when running the executable
             "alphafold" : \
                {# The settings for the output data frame
                 "df" : \
                    {# The settings for sorting the data frame
                     "sort" : \
                        {# The sorting criteria
                         "sort_by" : list,
                         # For each sorting criterion, whether
                         # to sortin ascending order
                         "ascending" : list},
                     # The settings for filtering the data frame
                     "filter" : dict}}},
        },

    #-----------------------------------------------------------------#

     # The template for the configuration for 'slimfast_get_homologs'
     "slimfast_get_homologs" : \

        {# The directory where the executable will be run
         "wd" : (str, os.getcwd()),
         
         # The settings to be used for the structures
         "structures" : \
            {# The settings to be used when running the executable
             "homolog" : \
                {# The settings for the output data frame
                 "df" : \
                    {# The settings for sorting the data frame
                     "sort" : \
                        {# The sorting criteria
                         "sort_by" : list,
                         # For each sorting criterion, whether
                         # to sortin ascending order
                         "ascending" : list},
                     # The settings for filtering the data frame
                     "filter" : dict}}},

         # The settings to run HHblits
         "hhblits" : \
            {# The options to run the program
             "run" : \
                {# The executable
                 "executable" : str,
                 # The name of the directory where to run the
                 # executable
                 "wd" : (str, ""),
                 # The name of the file where the output of the
                 # run will be stored
                 "outfile" : str,
                 # The options to launch the executable
                 "options" : dict}},

         # The settings to run HHsearch
         "hhsearch" : \
            {# The options to run the program
             "run" : \
                {# The executable
                 "executable" : str,
                 # The name of the directory where to run the
                 # executable
                 "wd" : (str, ""),
                 # The name of the file where the output of the
                 # run will be stored
                 "outfile" : str,
                 # The options to launch the executable
                 "options" : dict}},
        },

    #-----------------------------------------------------------------#

     # The template for the configuration for 'slimfast_get_homologs'
     "slimfast_find_slims" : \

        {# The directory where the executable will be run
         "wd" :  (str, os.getcwd()),
         
         # The settings to be used for the structures
         "structures" : \
            {# The minimum length (in number of residues) that a
             # disordered region needs to be to be considered a true
             # disordered region
             "min_len_disorder" : int,
             # # The name of the directory where to save the structures
             "struct_dir" : (str, "")},

         # The settings to run DSSP
         "dssp" : \
            {# The options to run the program
             "run" : \
                {# The executable
                 "executable" : str,
                 # The name of the directory where to run the
                 # executable
                 "wd" : (str, ""),
                 # The options to launch the executable
                 "options" : (dict, {})},
             # The processing options for the DSSP output
             "processing" : \
                {# The mapping of DSSP symbols to custom symbols
                 "dssp_mapping" : dict,
                 # How to represent unstructured residues
                 "whitespace_repr" : str,
                 # Which symbols of the DSSP mapping represent
                 # disordered residues
                 "disorder_symbols" : list}},

         # The settings to run NACCESS
         "naccess" : \
            {# The options to run the program
             "run" : \
                {# The executable
                 "executable" : str,
                 # The name of the directory where to run the
                 # executable
                 "wd" : (str, ""),
                 # The options to launch the executable
                 "options" : (dict, {})}},
        },

    }


#######################################################################


# The default level for logging
LOG_LEVEL = log.WARNING

# The options for logging
LOG_OPTIONS = \
   {# The format for log records
    "format" : "{asctime}:{levelname}:{name}:{message}",
    # The format for data and time in log records
    "datefmt" : "%Y-%m-%d,%H:%M",
    # The style of the log records
    "style" : "{"}

# The default names for the executables' log files
LOG_FILES = \
    {# The log file for 'slimfast_get_structures'
     "slimfast_get_structures" : \
        "slimfast_get_structures.log",

     # The log file for 'slimfast_get_alphafold'
     "slimfast_get_alphafold" : \
        "slimfast_get_alphafold.log",

     # The log file for 'slimfast_get_homologs'
     "slimfast_get_homologs" : \
        "slimfast_get_homologs.log",

     # The log file for 'slimfast_find_slims'
     "slimfast_find_slims" : \
        "slimfast_find_slims.log",

     # The log file for 'slimfast_run_pipeline'
     "slimfast_run_pipeline" : \
        "slimfast_run_pipeline.log"}


# The list of Dask loggers
DASK_LOGGERS = [\
    "distributed.batched",
    "distributed.client",
    "distributed.comm",
    "distributed.core",
    "distributed.diskutils",
    "distributed.http.proxy",
    "distributed.nanny",
    "distributed.scheduler",
    "distributed.worker",
    "bokeh.server"]

# The configuration for named formatters used when logging
CONFIG_FORMATTERS = \
    {# A generic formatter for log messages
     "generic_formatter" : \
    
        {# A format string in the given 'style' for the
         # logged output as a whole 
         "fmt" : LOG_OPTIONS["format"],

         # A format string in the given 'style' for the
         # date/time portion of the logged output
         "datefmt" : LOG_OPTIONS["datefmt"],

         # Can be one of '%', '{' or '$' and determines
         # how the format string will be merged with
         # its data
         "style" : LOG_OPTIONS["style"]},
    }

# The configuration for named filters used when logging
CONFIG_FILTERS = \
    {# Ignore messages starting with specific strings
     # from the 'distributed.worker' logger
     "distributed_filter" : \

        {# The custom 'factory' for the filter
         "()" : "SLiMfast.util.LevelContentFilter",
         
         # The level of the records that should be
         # ignored
         "level" : "INFO",

         # How the record should start to be ignored
         "start" : ["Starting Worker plugin",
                    "---",
                    "Local Directory:",
                    "Memory:",
                    "Threads:",
                    "Registered to:",
                    "Stopping worker",
                    "Waiting to connect to:",
                    "Listening to:",
                    "Start worker at:",
                    "Worker name:",
                    "dashboard at:",
                    "Worker closed",
                    "Found state lock",
                    "Starting established connection",
                    "Keep-alive"]},
    }


#######################################################################


# The default number of processes (cores) to be used when running
# 'slimfast_run_pipeline'
NPROC_DEFAULT = 1


#######################################################################


# The directory where the configuration files for running
CONFIG_RUN_DIR = \
        os.path.join(os.path.dirname(__file__), "SLiMfast", "config_run")

#######################################################################


# The columns of the NACCESS output whose data should be reported
NACCESS_COLS = ["abs_total_acc", "rel_total_acc"]


#######################################################################


# Name of the input file containing the SLiMs of interest for a
# given protein (it is automatically created by 
# 'slimfast_run_pipeline' from the input CSV file) and options
# to write it
INPUT_SLIMS_FILE = \
   {# The name of the file
    "name" : "slims.csv",
    # The options to write it
    "options" : \
      {# The field separator to be used
       "sep" : ",",
       # Whether to write out the index
       "index" : False,
       # Whether to write out the header
       "header" : True}}

# The default names and options for the output files produced by the
# SLiMfast executables
OUTPUT_FILES = \
    {# The CSV file produced by 'slimfast_get_structures'
     "protein" : \
        {# The name of the file
         "name" : "structures.csv",
         # The options used to write the file
         "options" : \
            {# The field separator
             "sep" : ",",
             # Whether to write out the index
             "index" : False,
             # The format to be used for dates
             "date_format" : "%Y/%m/%d"}},
     
     # The CSV file produced by 'slimfast_get_alphafold'
     "alphafold" : \
        {# The name of the file
         "name" : "alphafold.csv",
         # The options used to write the file
         "options" : \
            {# The field separator
             "sep" : ",",
             # Whether to write out the index
             "index" : False,
             # The format to be used for dates
             "date_format" : "%Y/%m/%d"}},

     # The CSV file produced by 'slimfast_get_homologs'
     "homologs" : \
        {# The name of the file
         "name" : "homologs.csv",
         # The options used to write the file
         "options" : \
            {# The field separator
             "sep" : ",",
             # Whether to write out the index
             "index" : False,
             # The format to be used for dates
             "date_format" : "%Y/%m/%d"}},

     # The JSON file produced by 'slimfast_get_homologs'
     "alignments" : \
        {# The name of the file
         "name" : "alignments.json"},

     # The CSV file produced by slimfast_find_slims'
     # containing the SLiMs found in the structures of
     # the protein of interest found in the Protein Data
     # Bank
     "slims_protein" : \
        {# The name of the file
         "name" : "slims_protein.csv",
         # The options used to write the file
         "options" : \
            {# The field separator
             "sep" : ",",
             # Whether to write out the index
             "index" : False}},

     # The CSV file produced by 'slimfast_find_slims'
     # containing the SLiMs found in the structures of
     # the protein of interest found in the AlphaFold
     # Protein Structure Database
     "slims_alphafold" : \
        {# The name of the file
         "name" : "slims_alphafold.csv",
         # The options used to write the file
         "options" : \
            {# The field separator
             "sep" : ",",
             # Whether to write out the index
             "index" : False}},

     # The output CSV file produced by 'slimfast_find_slims'
     # containing the SLiMs found in the structures of the
     # homologs of the protein of interest found in the Protein
     # Data Bank
     "slims_homologs" : \
        {# The name of the file
         "name" : "slims_homologs.csv",
         # The options used to write the file
         "options" : \
            {# The field separator
             "sep" : ",",
             # Whether to write out the index
             "index" : False}},

     # The output CSV file produced by 'slimfast_find_slims'
     # containing the SLiMs annotated through DisProt
     "slims_disprot" : \
        {# The name of the file
         "name" : "slims_disprot.csv",
         # The options used to write the file
         "options" : \
            {# The field separator
             "sep" : ",",
             # Whether to write out the index
             "index" : False}},
    }


#######################################################################


# The names of the different fields in the output data frames
DF_COLUMNS = \
    {# The name of the field containing a structure's ID
     "struct_id" : "struct_id",

     # The name of the field containing the ID of a chain
     # present in a structure
     "struct_chain" : "struct_chain",

     # The name of the field containing a structure's AlphaFold ID
     "alphafold_id" : "alphafold_id",

     # The name of the field containing the fragment number
     # for a structure coming from the AlphaFold Protein
     # Structure Database
     "struct_fragment" : "struct_fragment",

     # The name of the field containing the version of a structure
     # deposited in the AlphaFold Protein Structure Database
     "struct_version" : "struct_version",

     # The name of the field containing ranges representing
     # segments of a protein covered by a structure using the
     # structure's own residue numbering
     "struct_segments" : "struct_segments",

     # The name of the field containing the UniProt ID of a
     # protein present in a structure
     "uniprot_id" : "uniprot_id",

     # The name of the field containing the isoform of a protein
     # present in a structure
     "uniprot_isoform" : "uniprot_isoform",

     # The name of the field containing ranges representing
     # segments of a protein covered by a structure using the
     # UniProt residue numbering
     "uniprot_segments" : "uniprot_segments",

     # The name of the field containing the method used to obtain
     # the structure
     "method" : "method",

     # The name of the field containing the resolution of the
     # structure
     "resolution" : "resolution",

     # The name of the field containing a structure's deposition
     # date
     "deposition_date" : "deposition_date",

     # The name of the field containing the unique hit identifier
     # for a homolog
     "hit_id" : "hit_id"}

# The sorted columns for each of the output data frames
DF_SORTED_COLUMNS = \
    {# The sorted columns in the data frame produced by
     # 'slimfast_get_structures' 
     "protein" : \
        [DF_COLUMNS["struct_id"],
         DF_COLUMNS["struct_chain"],
         DF_COLUMNS["struct_segments"],
         DF_COLUMNS["uniprot_id"],
         DF_COLUMNS["uniprot_segments"],
         DF_COLUMNS["method"],
         DF_COLUMNS["resolution"],
         DF_COLUMNS["deposition_date"]],

     # The sorted columns in the data frame produced by
     # 'slimfast_get_alphafold'
     "alphafold" : \
        [DF_COLUMNS["struct_id"],
         DF_COLUMNS["struct_chain"],
         DF_COLUMNS["struct_segments"],
         DF_COLUMNS["struct_version"],
         DF_COLUMNS["struct_fragment"],
         DF_COLUMNS["uniprot_id"],
         DF_COLUMNS["uniprot_segments"],
         DF_COLUMNS["uniprot_isoform"],
         DF_COLUMNS["alphafold_id"],
         DF_COLUMNS["deposition_date"]],
    }


# The field separators to be used in the output data frames
# for different elements
DF_SEPARATORS = \
    {# The separator for numeric ranges
     "range" : "-",

     # The separator delimiting the starting and
     # ending points of different UniProt segments
     "segments" : ";"}
