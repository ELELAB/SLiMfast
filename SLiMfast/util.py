#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    util.py
#
#    Miscellanea utility functions for SLiMfast.
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
__doc__ = "Miscellanea utility functions for SLiMfast."


#######################################################################


# Standard library
import logging as log
import math
import os
import re
import requests
import sys
# Third-party packages
import pandas as pd
import yaml
# SLiMfast
from . import (
    alphafold_tools,
    config_tools,
    dataframe_tools,
    defaults,
    disprot_tools,
    dssp_tools,
    hhsearch_tools,
    naccess_tools,
    pdb_tools,
    uniprot_tools,
    run_tools)


#######################################################################


# Get the module's logger
logger = log.getLogger(__name__)


############################## CLASSES ################################


class LevelContentFilter(log.Filter):

    """
    A custom ``logging.Filter`` class to filter log records based
    on their level and content.
    """


    def __init__(self,
                 level = "WARNING",
                 start = None,
                 end = None,
                 content = None):
        """Initialize the filter.

        Parameters
        ----------
        level : ``str`` or ``logging`` level, ``logging.WARNING``
            The level at which and below which messages should
            be filtered out.

        start : ``list``, optional
            A list of strings. Log records starting with any
            of these strings will be filtered out.

            Note that the log record is stripped of any leading
            and traling blank spaces before checking whether it
            starts with any of the provided strings.

        end : ``list``, optional
            A list of strings. Log records ending with any
            of these strings will be filtered out.

            Note that the log record is stripped of any leading
            and trailing blank spaces before checking whether it
            ends with any of the provided strings.

        content : ``list``, optional
            A list of strings. Log records containing any
            of these strings will be filtered out.
        """

        # Call the parent class' initialization method
        super().__init__()

        # Save the level below which messages from the
        # specified loggers should be ignored
        self.level = log._nameToLevel[level]

        # Save the strings the log record needs to start with
        # to be ignored - if the log record starts with any
        # of the specified strings, it will be ignored
        self.start = start

        # Save the strings the log record needs to end with
        # to be ignored - if the log record ends with any
        # of the specified strings, it will be ignored
        self.end = end

        # Save the strings the log record needs to contain
        # to be ignored - if the log record contains any of
        # the specified strings, it will be ignored
        self.content = content


    def filter(self,
               record):
        """Decide whether a log record will be emitted or not.

        Parameters
        ----------
        record : ``logging.LogRecord``
            The log record.

        Returns
        -------
        emit_record : ``bool``
            Whether the log record should be emitted or not.
        """

        # If the record's level is below or equal to the selected
        # level
        if record.levelno <= self.level:

            # If we specified strings the record needs to start with
            # to be ignored
            if self.start is not None:

                # For each string
                for string in self.start:

                    # If the record starts with the selected string
                    if record.msg.strip().startswith(string):

                        # Do not log the record
                        return False

                # If we reached this point, log the record
                return True

            #---------------------------------------------------------#

            # If we specified strings the record needs to end with
            # to be ignored
            if self.end is not None:

                # For each string
                for string in self.end:

                    # If the record ends with the selected string
                    if record.msg.strip().endswith(string):

                        # Do not log the record
                        return False

                # If we reached this point, log the record
                return True

            #---------------------------------------------------------#

            # If we specified strings the record needs to contain
            # to be ignored
            if self.content is not None:

                # For each string
                for string in self.content:

                    # If the record contains the selected string
                    if string in record:

                        # Do not log the record
                        return False

                # If we reached this point, log the record
                return True

        #-------------------------------------------------------------#

        # If we reached this point, log the record
        return True


########################## GENERIC UTILITIES ##########################


def get_wd(wd):
    """Get the working directory where a step of the pipeline or
    the whole pipeline will be run.

    Parameters
    ----------
    wd : ``str``
        The working directory as specified in the configuration file.

    Returns
    -------
    wd : ``str``
        The path to the working directory. 
    """

    # If the working directory is defined as an empty string
    if wd is None or wd == "":
        
        # The working directory will be the current working directory
        wd = os.getcwd()

    #-----------------------------------------------------------------#
    
    # Otherwise
    else:
        
        # If it is the name of a directory
        if os.path.basename(wd) == wd:
            
            # The working directory will be a new directory with than
            # name created inside the current working directory
            wd = os.path.join(os.getcwd(), wd)

            # Create the working directory
            os.mkdir(wd)

            # Inform the user that the directory was created
            infostr = \
                f"The directory '{os.path.basename(wd)}' was " \
                f"created inside '{os.getcwd()}'."
            log.info(infostr)
        
        # If it is the path to an existing directory
        else:
            
            # Get the absolute path of the directory and set
            # it as the working directory - assume that the
            # directory already exists
            wd = os.path.abspath(wd)

    #-----------------------------------------------------------------#

    # Return the path to the working directory
    return wd


def get_abspath(path = None):
    """Given a path, return its absolute path.

    Return ``None`` if the given path is ``None``.

    Parameters
    ----------
    path : ``str``, optional
        The original path.

    Returns
    -------
    abs_path : ``str`` or ``None``
        The absolute path of the original path if ``path`` was
        not ``None``; otherwise, ``None``. 
    """

    # Return the absolute path or None
    return os.path.abspath(path) if path is not None else path


def get_values_from_dict(d,
                         keys,
                         expected_num_values = None,
                         insert_first = None):
    """Get a list of values from a dictionary, given their keys.
    
    Parameters
    ----------
    d : ``dict``
        The input dictionary.

    keys : ``list``
        A list of keys whose associated values need to be retrieved.

    expected_num_values : ``int``, optional
        The expected number of values for each key.

        If not provided, one value per key is expected.

        If the key does not exist, ``None`` will be added to the
        final list for as many times as the expected number of values.

    insert_first : ``list``, optional
        A list of elements to be inserted first, if any.

    Returns
    -------
    values_list : ``list``
        The list of values. 
    """

    # Create an empty list of itens
    items = []

    #-----------------------------------------------------------------#

    # If there are elements to be inserted first
    if insert_first is not None:

        # Add them
        items.extend(insert_first)

    #-----------------------------------------------------------------#

    # If there is no expected number of values
    if expected_num_values is None:

        # Set it to 1 for each key
        expected_num_values = [1 for k in keys]

    #-----------------------------------------------------------------#

    # For each key
    for k, num_v in zip(keys, expected_num_values):

        # If the key is in the dictionary
        if k in d.keys():
            
            # Get the corresponding value
            v = d[k]

            # If the value is a list or a tuple
            if isinstance(v, (list, tuple)):

                # Extend the list
                items.extend(v)

            # Otherwise
            else:

                # Append the value to the list
                items.append(v)

        # Otherwise
        else:

            # Append None to the list for as many times as the
            # expected number of values
            items.extend([None] * num_v)

    #-----------------------------------------------------------------#

    # Return the list
    return items


########################## LOGGING UTILITIES ##########################


def get_log_handlers(log_to_stdout = True,
                     log_file_handler_class = None,
                     log_file_handler_options = None,
                     log_level = defaults.LOG_LEVEL):
    """Get the handlers to set up the logging.

    Parameters
    ----------
    log_to_stdout : ``bool``, ``True``
        Whether to write log messages to the standard output.

    log_file_handler_class : ``logging.FileHandler``, optional
        A ``FileHandler`` class to construct the handler that
        will log to a log file.

        If not provided, the log messages will not be written
        to any file.

    log_file_handler_options : ``dict``, optional
        A dictionary of options to set up the handler that will
        log to a file.

        It must be provided if ``log_file_class`` is provided.

    log_level : ``str``, optional
        The level below which log messages will be silenced.

        By default, it is the one defined by 
        ``SLiMfast.defaults.LOG_LEVEL``.
    """

    # Create a list to store the handlers
    handlers = []

    #-----------------------------------------------------------------#

    # If the user requested logging to the standard output
    if log_to_stdout:

        # Create the StreamHandler
        handler = log.StreamHandler()

        # Set the handler's level
        handler.setLevel(log_level)

        # Add the handler to the list of handlers
        handlers.append(handler)

    #-----------------------------------------------------------------#

    # If the user specified a class to create a log file
    if log_file_handler_class is not None:

        # Try to get the path to the log file
        try:

            log_file = log_file_handler_options["filename"]

        # If the user did not pass any options
        except TypeError as e:

            # Raise an error
            errstr = \
                "If 'log_file_class' is provided, a dictionary of " \
                "'log_file_options' must be passed as well."
            raise TypeError(errstr)

        # If the user did not pass the 'filename' option
        except KeyError as e:

            # Raise an error
            errstr = \
                "A 'filename' must be provided in the dictionary of " \
                "'log_file_options."
            raise KeyError(errstr)

        # If the file already exists
        if os.path.exists(log_file):

            # Remove it
            os.remove(log_file)

        #-------------------------------------------------------------#

        # Create the corresponding handler
        handler = log_file_handler_class(**log_file_handler_options)

        # Set the handler's level
        handler.setLevel(log_level)

        # Add the handler to the list of handlers
        handlers.append(handler)

    #-----------------------------------------------------------------#
    
    # Return the handlers
    return handlers


def get_dask_logging_config(log_to_stdout = True,
                            log_file = None,
                            log_level = "ERROR"):
    """Get the logging configuration for Dask/distributes
    loggers.

    Parameters
    ----------
    log_to_stdout : ``bool``, ``True``
        Whether to log messages to the standard output.

    log_file : ``str``, optional
        The name of the log file where to write the
        log messages.

        If not provided, log messages will not be
        written to any file.

    log_level : ``str``, ``"ERROR"``
        Below which level log messages should be silenced.

    Returns
    -------
    dask_logging_config : ``dict``
        The logging configuration for Dask/distributed
        loggers.
    """

    # Initialize the Logging configuration - it follows the
    # "configuration dictionar schema" provided in
    # https://docs.python.org/3/library/logging.config.html
    # #configuration-dictionary-schema
    dask_logging_config = {\
        
        # Version of the logging configuration - so far,
        # only 1 is supported
        "version": 1,
        
        # The formatters for the log records - each key
        # represents the name of a 'logging.Formatter' object
        # and the dictionary associated with it contains the
        # options to initialize it
        "formatters" : {\
            
            # Use the generic formatter
            "generic_formatter" : \
                defaults.CONFIG_FORMATTERS["generic_formatter"],
        }, 

        # The filters to be used for log records - each key
        # represents the name of a 'logging.Filter' object
        # and the dictionary associated with it contains the
        # options to initialize it
        "filters" : {\
            
            # Use the custom filter for the 'distributed.worker'
            # logger
            "distributed_filter" : \
                defaults.CONFIG_FILTERS["distributed_filter"],
        },
       
        # The handlers to be used when logging - each key
        # represents the name of a 'logging.Handler' object
        # and the dictionary associated with it contains the
        # options to initialize it
        "handlers": {
            
            # A handler to log to a rotating file
            "rotating_file_handler": \
                {# The class of the handler to initialize
                 "class": "logging.handlers.RotatingFileHandler",
                 # 'filename' needs to be set before using the
                 # configuration
                 "filename": log_file,
                 # The formatter for the log records
                 "formatter" : "generic_formatter"},
            
            # A handler to log to the standard output
            "stream_handler": \
                {# The class of the handler to initialize
                 "class": "logging.StreamHandler",
                 # The formatter for the log records
                 "formatter" : "generic_formatter"},
        },

        # The loggers to configure - each key represents
        # the name of a 'logging.Logger' object and the
        # dictionary associated with it contains the
        # options to configure it. The dictionary of
        # loggers is empty because we are going to fill
        # it later
        "loggers" : {},
    }

    #-----------------------------------------------------------------#

    # Initnialize an empty list to store the handlers to be used
    handlers = []

    # If the user requested logging to the standard output
    if log_to_stdout:

        # Add the corresponding handler to the list
        handlers.append("stream_handler")

    # Otherwise
    else:

        # Remove the corresponding handler from the configuration
        dask_logging_config["handlers"].pop("stream_handler")

    # If the user requested logging to a file
    if log_file is not None:

        # Add the corresponding handler to the list
        handlers.append("rotating_file_handler")

    # Otherwise
    else:

        # Remove the corresponding handler from the configuration
        dask_logging_config["handlers"].pop("rotating_file_handler")

    #-----------------------------------------------------------------#

    # For each Dask logger
    for logger in defaults.DASK_LOGGERS:
  
        # Configure it
        dask_logging_config["loggers"][logger] = \
            {"level" : log_level,
             "filters" : ["distributed_filter"],
             "handlers" : handlers}

    #-----------------------------------------------------------------#

    # Return the configuration
    return dask_logging_config


####################### FILE-LOADING UTILITIES ########################


def load_config(config_file,
                config_type):
    """Load the configuration from a YAML file.

    Parameters
    ----------
    config_file : ``str``
        The YAML configuration file.

    config_type : ``str``, {``"slimfast_get_structures"``,
        ``"slimfast_get_alphafold"``, ``"slimfast_get_homologs"``,
        ``"slimfast_find_slims"``, ``"slimfast_run_pipeline"``}
        The type of configuration we are loading.

    Returns
    -------
    config : ``dict``
        The configuration.
    """

    # Get the name of the configuration file without extension
    # and the extension
    config_name, config_ext = \
        os.path.splitext(os.path.basename(config_file))

    #-----------------------------------------------------------------#

    # If the configuration file is a name without extension
    if config_file == config_name:

        # Inform the user that we are looking for the configuration
        # file in the default directory
        infostr = \
            f"Looking for '{config_name}' in " \
            f"'{defaults.CONFIG_RUN_DIR}'..."
        logger.info(infostr)

        # Assume it is a file with that name inside the default
        # directory for configuration files
        config_file = \
            os.path.join(defaults.CONFIG_RUN_DIR,
                         config_name + config_ext)

    #-----------------------------------------------------------------#

    # Otherwise
    else:

        # Get the absolute path to the configuration file
        config_file = os.path.abspath(config_file)

    #-----------------------------------------------------------------#

    # Load the configuration
    config = yaml.safe_load(open(config_file, "r"))

    #-----------------------------------------------------------------#

    # Get the working directory where the step/pipeline will be run
    wd = get_wd(wd = config.get("wd"))

    # Set it
    config["wd"] = wd

    #-----------------------------------------------------------------#

    # Check the configuration
    config_tools.check_against_template(config = config,
                                        config_type = config_type)

    #-----------------------------------------------------------------#

    # Return the configuration
    return config


def load_slims(slims_file):
    """Load the SLiMs of interest from a CSV file and return a list
    containing them.

    The file is a expected to contain a comma-separated data frame
    with the following mandatory fields:

    - ``slim_name`` : a string representing the SLiM's name.
    - ``slim_start`` : an integer representing where the SLiM
      starts in the protein sequence.
    - ``slim_seq`` : a string representing the SLiM's sequence.
    - ``slim_end`` : an integer representing where the SLiM ends
      in the protein sequence.

    Additionally, one more field may be specified:

    - ``slim_acc_res`` period-separated list of integers representing
      the positions of the residues in the SLiM on which the analysis
      of solvent accessibility will be performed. The positions should
      be specified relative to the full UniProt sequence of the protein
      the SLiM belongs to.

    Parameters
    ----------
    slims_file : ``str``
        The file containing the SLiMs.

    Returns
    -------
    slims : ``list``
        A list containing the SLiMs found in the file. Each SLiM is
        represented by a list containing:

        - A string representing the SLiM's name.
        - An integer representing where the SLiM starts in the protein
          sequence.
        - A string representing the SLiM's sequence.
        - An integer representing where the SLiM ends in the protein
          sequence.
        - If a list of residues on which solvent accessibility should
          be calculated was passed, a tuple containing two tuples,
          the first containing a list of integers representing the
          absolute sequence positions of the residues, and the second
          containing a list of integers representing the relative
          sequence positions of the residues (with respect to the
          SLiM's sequence). If no list of residues was passed, the
          tuple will contain two ``None`` elements. 
    """

    # Set the list of columns that can be found in the file
    VALID_COLS = \
        ["slim_name", "slim_start", "slim_seq", "slim_end",
         "slim_acc_res"]

    # Get a string representing the valid columns
    valid_cols_str = ", ".join([f"'{c}'" for c in VALID_COLS])

    #-----------------------------------------------------------------#

    # Load the SLiMs as a data frame
    df_slims = pd.read_csv(slims_file,
                           sep = ",")

    #-----------------------------------------------------------------#

    # Convert the data frame to a list of tuples
    slims_tuples = list(df_slims.itertuples(index = False,
                                            name = None))
    #-----------------------------------------------------------------#

    # Create an empty list to store the results
    slims = []

    #-----------------------------------------------------------------#

    # If there are four columns
    if len(slims_tuples[0]) == 4:

        # If they are not the supported ones
        if df_slims.columns.tolist() != VALID_COLS[:-1]:

            # Set a string representing the invalid columns
            cols_str = \
                ", ".join([f"'{c}'" for c in df_slims.columns])

            # Raise an error
            errstr = \
                "Invalid columns were found in " \
                f"'{slims_file}': {cols_str}. Valid columns " \
                f"are: {valid_cols_str}."
            raise ValueError(errstr)

        # For each tuple
        for tup in slims_tuples:

            # Add the tuple to the final list, plus one field
            # to indicate the (missing) residues for which to
            # calculate the solvent accessibility
            slims.append(\
                (tup[0], tup[1], tup[2], tup[3], (None, None)))

    #-----------------------------------------------------------------#

    # If there are five columns
    elif len(slims_tuples[0]) == 5:

        # If they are not the supported ones
        if df_slims.columns.tolist() != VALID_COLS:

            # Set a string representing the invalid columns
            cols_str = \
                ", ".join([f"'{c}'" for c in df_slims.columns])

            # Raise an error
            errstr = \
                "Invalid columns were found in " \
                f"'{slims_file}': {cols_str}. Invalid columns " \
                f"are: {valid_cols_str}."
            raise ValueError(errstr)

        # For each tuple
        for tup in slims_tuples:

            # Get the start of the SLiM
            slim_start = tup[1]

            # Get the residues for which we have to calculate the
            # solvent accessibility
            slim_acc_res_abs = \
                [int(item) for item in tup[4].split(":")]

            # Get the relative positions of the residues with respect
            # to the SLiM
            slim_acc_res_rel = \
                [res - slim_start for res in slim_acc_res_abs]

            # Add a new tuple with the updated data to the final list
            slims.append(\
                (tup[0], tup[1], tup[2], tup[3], 
                 (tuple(slim_acc_res_abs), tuple(slim_acc_res_rel))))

    #-----------------------------------------------------------------#

    # Return the list of SLiMs
    return slims


def load_ids_slims(ids_slims_file):
    """Load the UniProt IDs of the proteins of interest and the
    associated SLiMs from a CSV file and return a dictionary mapping
    each protein's UniProt ID to a data frame containing the SLiMs of
    interest for that protein.

    Parameters
    ----------
    ids_slims_file : ``str``
        The input CSV file.

    Returns
    -------
    id2slims : ``dict``
        A dictionary mapping each protein's UniProt ID to a data
        frame containing the SLiMs of interest for that protein.
    """

    # Set the field that should contain the UniProt IDs of the
    # proteins of interest
    uniprot_id_col = defaults.DF_COLUMNS["uniprot_id"]

    #-----------------------------------------------------------------#

    # Read in the CSV file as a data frame
    df = pd.read_csv(ids_slims_file,
                     sep = ",")

    #-----------------------------------------------------------------#

    # If there are less than the expected five columns in the data
    # frame
    if len(df.columns) < 5:

        # If there is only one column and it contains the UniProt
        # IDs of the proteins of interest (which can happen when
        # the SLiM search is not requested as part of the SLiMfast
        # pipeline)
        if len(df.columns) == 1 \
        and uniprot_id_col in df.columns.tolist():

            # Return a data frame mapping the UniProt IDs to None
            return {k : None for k in df[uniprot_id_col]}

        # Otherwise
        else:

            # Raise an error
            errstr = \
                f"The file '{ids_slims_file}' contains fewer than " \
                "the expected five fields. The only case in which " \
                f"this is acceptable is when the data frame only " \
                f"contains the column '{uniprot_id_col}' with the " \
                "UniProt IDs of the proteins of interest."
            raise ValueError(errstr)

    #-----------------------------------------------------------------#

    # Convert the data frame into a dictionary keyed on the unique
    # UniProt IDs (SLiMs for the same IDs are, therefore, grouped
    d = df.groupby(uniprot_id_col).agg(\
            lambda x: x.tolist()).to_dict("index")

    #-----------------------------------------------------------------#

    # Return a dictionary where keys are the UniProt IDs and values
    # are data frames containing the SLiMs for each UniProt ID
    return {k : pd.DataFrame(v) for k, v in d.items()}


def get_structures_df(uniprot_id):
    """Given a protein's UniProt ID, return a data frame containing
    information about the experimental structures of the protein
    available in the Protein Data Bank.

    Parameters
    ----------
    uniprot_id : ``str``
        The protein's UniProt ID.

    Returns
    -------
    structures_df : ``pandas.DataFrame``
        A data frame containing information about the available
        experimental structures found in the Protein Data Bank for
        the protein.  
    """

    # Set the separator for ranges
    range_sep = defaults.DF_SEPARATORS["range"]

    # Set the separator for protein segments
    segments_sep = defaults.DF_SEPARATORS["segments"]

    #-----------------------------------------------------------------#

    # Get the name of the column that will contain the structures' IDs
    struct_id_col = defaults.DF_COLUMNS["struct_id"]

    # Get the name of the column that will contain the structures'
    # chains
    struct_chain_col = defaults.DF_COLUMNS["struct_chain"]

    # Get the name of the column that will contain the protein's
    # segments covered by the structures, in the structures' own
    # numbering
    struct_segments_col = defaults.DF_COLUMNS["struct_segments"]

    # Get the name of the column that will contain the protein's
    # UniProt ID
    uniprot_id_col = defaults.DF_COLUMNS["uniprot_id"]

    # Get the name of the column that will contain the protein's
    # segments covered by the structures, in UniProt numbering
    uniprot_segments_col = defaults.DF_COLUMNS["uniprot_segments"]

    #-----------------------------------------------------------------#
    
    # Create an empty list where to store the data (it will be
    # converted into a data frame)
    data = []

    #-----------------------------------------------------------------#

    # Get the IDs of the available experimental structures for the
    # protein of interest in the Protein Data Bank (it yields one
    # structure at a time)
    pdbs = uniprot_tools.get_pdbs(uniprot_id = uniprot_id)


    #-----------------------------------------------------------------#
    
    # For each ID
    for struct_id in pdbs:
        
        # Get the mappings between the segments of the protein found
        # in the PDB structure and the segments of the full protein
        # sequence they correspond to
        mappings = \
            pdb_tools.get_protein_segments(struct_id = struct_id,
                                           uniprot_id = uniprot_id)

        #-------------------------------------------------------------#
        
        # For each chain in the PDB structure and the corresponding
        # segments of the full UniProt sequence covered by that chain
        for struct_chain, segments in mappings.items():
            
            # Get the starting and ending points of the segments of
            # the protein that are covered by the current chain in
            # the PDB structure, using the UniProt numbering
            unp_segments = \
                [range_sep.join(\
                    [str(t) for t in s[1]]) for s in segments]

            # Join the segments using the UniProt numbering
            # in a formatted string
            unp_segments_joined = segments_sep.join(unp_segments)

            #---------------------------------------------------------#

            # Get the starting and ending points of the segments of
            # the protein that are covered by the current chain in
            # the PDB structure, using the numbering found in the
            # PDB structure
            struct_segments = \
                [range_sep.join(\
                    [str(t) for t in s[0]]) for s in segments]

            # Join the segments using the PDB numbering in a
            # formatted string
            struct_segments_joined = segments_sep.join(struct_segments)

            #---------------------------------------------------------#
            
            # Add the data of interest to the final list
            data.append(\
                {struct_id_col : struct_id,
                 struct_chain_col : struct_chain,
                 struct_segments_col : struct_segments_joined,
                 uniprot_id_col : uniprot_id,
                 uniprot_segments_col : unp_segments_joined})

    #-----------------------------------------------------------------#

    # If no matching PDB structures were found for the UniProt ID
    if not data:
        
        # Return an empty data frame
        return pd.DataFrame()

    #-----------------------------------------------------------------#

    # Create a data frame from the data
    df = pd.DataFrame(data)

    #-----------------------------------------------------------------#

    # Return the data frame
    return df


def get_alphafold_df(uniprot_id):
    """Given a protein's UniProt ID, return a data frame containing
    information about the structures of the protein available in the
    AlphaFold Protein Structure Database.

    Parameters
    ----------
    uniprot_id : ``str``
        The protein's UniProt ID.

    Returns
    -------
    alphafold_df : ``pandas.DataFrame``
        A data frame containing information about the available
        structures found in the AlphaFold Protein Structure Database
        for the protein.
    """

    # Get the name of the column that will contain the structures' IDs
    struct_id_col = defaults.DF_COLUMNS["struct_id"]

    # Get the name of the column that will contain the protein's
    # UniProt ID
    uniprot_id_col = defaults.DF_COLUMNS["uniprot_id"]

    # Get the name of the column that will contain the protein's
    # fragment covered by each structures (for structures split over
    # different fragments)
    strcut_fragment_col = defaults.DF_COLUMNS["struct_fragment"]

    # Get the name of the column that will contain each structure's
    # version
    struct_version_col = defaults.DF_COLUMNS["struct_version"]

    #-----------------------------------------------------------------#
    
    # Create an empty list where to store the data (it will be
    # converted into a data frame)
    data = []

    #-----------------------------------------------------------------#

    # Get the IDs of the structures of the protein available in the
    # AlphaFold Protein Structure Database
    structs = \
        alphafold_tools.get_alphafold_data(uniprot_id = uniprot_id)

    #-----------------------------------------------------------------#

    # If no structures were found
    if not structs:
        
        # Return an empty data frame
        return pd.DataFrame()

    #-----------------------------------------------------------------#
    
    # For each structure found
    for struct_id, (version, fragment, _) in structs.items():
            
        # Append the data of interest to the dictionary
        data.append({struct_id_col : struct_id,
                     uniprot_id_col : uniprot_id,
                     strcut_fragment_col : fragment,
                     struct_version_col : version})

    #-----------------------------------------------------------------#

    # Create a data frame from the data
    df = pd.DataFrame(data)

    #-----------------------------------------------------------------#

    # Return the data frame
    return df


def add_structures_metadata(df,
                            source):
    """Given a data frame containing data about protein structures,
    add metadata about the structures to the data frame.

    Parameters
    ----------
    df : ``pandas.DataFrame``
        The data frame containing data about the structures.

    source : ``str``, {``"pdb"``, ``"alphafold"``}
        The source of the structures.

    Returns
    -------
    updated_df : ``pandas.DataFrame``
        An updated version of the data frame containing also the
        structure's metadata.
    """

    # Get the name of the column that contains the structures' IDs
    struct_id_col = defaults.DF_COLUMNS["struct_id"]

    #-----------------------------------------------------------------#

    # Create an empty list to store the metadata found for each
    # structure
    metadata_list = []

    #-----------------------------------------------------------------#

    # For each unique ID
    for struct_id in df[struct_id_col].unique():

        # If the structure comes from the Protein Data Bank
        if source == "pdb":

            # Try to get the structure's metadata
            try:
                
                metadata = \
                    pdb_tools.get_structure_metadata(\
                        struct_id = struct_id)

            # If something went wrong       
            except Exception as e:

                # Log a warning
                warnstr = \
                    "It was not possible to retrieve the " \
                    "metadata for the structure with ID " \
                    f"'{struct_id}' from the Protein Data " \
                    f"Bank. This structure will be ignored. " \
                    f"Error: {e}"
                logger.warning(warnstr)
                
                # Go to the next structure
                continue

        #-------------------------------------------------------------#

        # If the structure comes from the AlphaFoldDB
        elif source == "alphafold":
        
            # Try to get the structure's metadata
            try:
                
                metadata = \
                    alphafold_tools.get_structure_metadata(\
                        struct_id = struct_id)

            # If something went wrong       
            except Exception as e:

                # Log a warning
                warnstr = \
                    "It was not possible to fetch the " \
                    "metadata for the structure with ID " \
                    f"'{struct_id}' from the AlphaFold " \
                    f"Protein Structure Database. This structure " \
                    f"will be ignored. Error: {e}"
                logger.warning(warnstr)
                
                # Go to the next structure
                continue

        #-------------------------------------------------------------#

        # Add a field for the structure ID (it is needed to merge
        # the dataframes later)
        metadata.update({struct_id_col : struct_id})

        # Add the dictionary of metadata to the list
        metadata_list.append(metadata)

        # Inform the user that the metadata were collected
        # successfully for the current structure
        infostr = \
            f"The metadata for the structure with ID '{struct_id}' " \
            "were successfully collected."
        logger.info(infostr)

    #-----------------------------------------------------------------#

    # If metadata for at least one of the structures were retrieved
    if metadata_list:

        # Create a data frame containing the metadata for all
        # structures
        metadata_df = pd.DataFrame(metadata_list)

        # Combine the concatenated data frame with the data frame
        # containing the structures' metadata
        updated_df = pd.merge(left = df,
                              right = metadata_df,
                              on = struct_id_col,
                              how = "outer")

        # Return the updated data frame
        return updated_df

    #-----------------------------------------------------------------#

    # Otherwise
    else:

        # Inform the user that no metadata were found for the
        # structures
        warnstr = \
            "No metadata could be retrieved for any structure."
        logger.warning(warnstr)

        # Return the original data frame
        return df


def region_in_protein_segments(seq_range,
                               uniprot_segments,
                               struct_segments):
    """Check if a sequence belongs to at least one UniProt segment,
    given the sequence's numeric range.

    Parameters
    ----------
    seq_range : ``tuple``
        A tuple of integers representing the range covered by the
        current sequence in the context of the full protein sequence.
        
        Both ends of the range are included.

    uniprot_segments : ``str``
        A string representing the portions of the protein sequence
        covered by a structure using the UniProt numbering.

    struct_segments : ``str``
        A string representing the portions of the protein sequence
        covered by a structure using the structure's own numbering.

    Returns
    -------
    is_in_segments : ``bool``
        Whether the specified ``seq_range`` falls in any of the
        ``uniprot_segments`` covered by the structure. 
    """

    # Set the separator for ranges
    range_sep = defaults.DF_SEPARATORS["range"]

    # Set the separator for protein segments
    segments_sep = defaults.DF_SEPARATORS["segments"]

    #-----------------------------------------------------------------#

    # Unpack the range
    start_s, end_s = seq_range

    #-----------------------------------------------------------------#

    # Try to get the numeric ranges of all UniProt segments passed
    try:
        
        uniprot_ranges = \
            [(int(s.split(range_sep)[0]), int(s.split(range_sep)[1])) \
             for s in uniprot_segments.split(segments_sep)]

    # If something went wrong
    except ValueError:

        # Set the error message
        errstr = \
            "It was not possible to find the segments convered by " \
            "the structure (expressed in UniProt numbering), since " \
            "at least one extreme of at least one of the intervals " \
            "defining the segments is None."

        # Raise a new exception
        raise ValueError(errstr)

    #-----------------------------------------------------------------#

    # Try to get the numeric ranges of all segments of the structure
    # passed
    try:
        
        struct_ranges = \
            [(int(s.split(range_sep)[0]), int(s.split(range_sep)[1])) \
             for s in struct_segments.split(segments_sep)]

    # If something went wrong
    except ValueError:

        # Log the error
        errstr = \
            "It was not possible to find the segments covered " \
            "by the structure (expressed in the structure's own " \
            "numbering) since at least one extreme of at least one " \
            "one of the intervals defining the segments is None."
        
        # Raise a new exception
        raise ValueError(errstr)

    #-----------------------------------------------------------------#

    # For each UniProt segment
    for uniprot_start, uniprot_end in uniprot_ranges:

        # If the SLiM is found within the segment
        if start_s >= uniprot_start and end_s <= uniprot_end:
            
            # Return True                   
            return True

    #-----------------------------------------------------------------#

    # If the sequence was never found, return False
    return False


def region_in_homolog_segments(pdb_data,
                               struct_id,
                               struct_chain,
                               hit_al,
                               seq_range):
    """Check if a region of the protein of interest is found in the
    segments corresponding to the sequences aligned to such a
    region in the proteins' homologs.

    Parameters
    ----------
    pdb_data : ``list``
        A list with the lines parsed from a PDB file.

    struct_id : ``str``
        The structure's ID.

    struct_chain : ``str``
        The structure's chain we are currently investigating.

    hit_al : ``dict``
        A dictionary containing an alignment parsed from the
        output of HHsearch.

    seq_range : ``list`` or ``tuple``
        The starting and ending point of the protein region to
        find in the alignment.

    Returns
    -------
    is_in_segments : ``bool``
        Whether the specified ``seq_range`` falls in any of the
        segments` of the protein of interest covered by the
        structure. 

    start_in_t : ``int``
        The start of the SLiM of interest in the homolog, if
        ``is_in_segments`` is ``True``; otherwise, ``None``.

    slim_in_t : ``seq``
        The sequence of the SLiM of interest in the homolog, if
        ``is_in_segments`` is ``True``; otherwise, ``None``.

    end_in_t : ``int``
        The end of the SLiM of interest in the homolog, if
        ``is_in_segments`` is ``True``; otherwise, ``None``.
    """

    # Try to find the sequence in the homolog aligned to the
    # sequence of interest in the HHsearch alignment
    try:
        
        start_in_t, slim_in_t, end_in_t, is_match = \
            hhsearch_tools.get_aligned_sequence(hit_al = hit_al,
                                                seq_range = seq_range)
    
    # If something went wrong
    except Exception as e:

        # Set the error message
        start, end = seq_range
        errstr = \
            "It was not possible to find the sequence " \
            f"{start}-{end} in the alignment provided."
        
        # Re-raise the exception
        raise Exception(errstr) from e

    #-----------------------------------------------------------------#

    # If the aligned sequence does not perfectly match the sequence
    # of interest
    if not is_match:

        # Return
        return False, None, None, None

    #-----------------------------------------------------------------#

    # Get the numeric range covered by the aligned sequence in the
    # structure of the homolog
    slim_range = set(list(range(start_in_t, end_in_t)))

    #-----------------------------------------------------------------#

    # Try to get the missing residues in the structure
    try:
        
        missing_res = \
            pdb_tools.get_missing_residues(pdb_data = pdb_data)
    
    # If something went wrong
    except Exception as e:

        # Set the error message
        errstr = \
            f"Could not get the missing residues for the PDB " \
            f"structure with ID {struct_id}. Exception: {e}"
        
        # Re-raise the exception
        raise Exception(errstr) from e

    #-----------------------------------------------------------------#
    
    # Get the missing residues for the chain of interest
    missing_res_chain = set(missing_res[struct_chain])

    # If none of the residues in the aligned sequence of interest
    # is missing
    if not set(slim_range).intersection(missing_res_chain):

        # Since none of the SLiM residues were missing, return that
        # the region was found in the structure, plus the start,
        # sequence, and end of the SLiM in the region
        return True, start_in_t, slim_in_t, end_in_t

    #-----------------------------------------------------------------#

    # Return None for each of the structure-related info if
    # some residues of the aligned sequence of interest were
    # missing in the structure of the homolog
    return True, None, None, None


def find_slims(source,
               structures_df,
               slims_to_find,
               struct_dir,
               als = None):
    """Find SLiMs of interest in a set of structures.

    Parameters
    ----------
    source : ``str``, {``"structures"``, ``"alphafold"``,
        ``"homologs"``}
        Where the structures come from.

    structures_df : ``pandas.DataFrame``
        A data frame containing information about the structures.

    slims_to_find : ``list``
        A list containing the SLiMs for which DisProt annotations
        should be retrieved. Each SLiM should be represented by a
        list containing:

        - A string representing the SLiM's name.
        - An integer representing where the SLiM starts in the protein
          sequence.
        - A string representing the SLiM's sequence.
        - An integer representing where the SLiM ends in the protein
          sequence.
        - If a list of residues on which solvent accessibility should
          be calculated was passed, a tuple containing two tuples,
          the first containing a list of integers representing the
          absolute sequence positions of the residues, and the second
          containing a list of integers representing the relative
          sequence positions of the residues (with respect to the
          SLiM's sequence). If no list of residues was passed, the
          tuple will contain two ``None`` elements. 

    struct_dir : ``str``
        The directory where to save the structures' PDB files.

    als : ``dict``, optional
        The alignments associated with the structures, if ``source``
        is ``"homologs"``.

    Returns
    -------
    structures_and_slims : ``list``
        A list of lists containing as many entries as the unique
        structures and associated unique chains in which SLiMs
        were found.

        Each entry consists of a list containing the following
        elements:

        - A string represented the type of 'source' represented
          by the structure. This can be ``"structures"``,
          ``"alphafold"``, or ``"homologs"`` depending on whether
          the structure is an experimental structure of the
          protein of interest from the Protein Data Bank, a
          structure of the protein of interest from the AlphaFold
          Protein Structure Database,  or the experimental structure
          of a homolog of the protein of interest from the Protein
          Data Bank.

        - A string representing the ID of the structure. This
          is the PDB ID for experimental structures coming from
          the Protein Data Bank and the AlphaFold ID for
          structures coming from the AlphaFold Protein Structure
          Database.

        - A string representing the chain ID of the structure's
          chain the SLiMs were found in.

        - A string representing the path to the PDB file containing
          the structure.

        - A list representing the SLiMs found in the given
          structure's chain.
    """
        
    # Create an empty dictionary storing the data for the structures
    # of interest
    struct_data_dict = {}

    #-----------------------------------------------------------------#

    # Create an empty list to store the PDB files and the SLiMs to be
    # checked for disorder
    structures_and_slims = []

    #-----------------------------------------------------------------#

    # Create a set to store the SLiMs found
    slims_found = set()

    #-----------------------------------------------------------------#

    # For each row in the data frame (= structure)
    for index, row in structures_df.iterrows():

        #-------------------------------------------------------------#

        # Create an empty list to store the SLiMs found in each
        # structure
        slims_in_struct = []

        #-------------------------------------------------------------#

        # Get the structure ID, chain, the segments of the protein
        # of interest covered by the structure
        struct_id = row[defaults.DF_COLUMNS["struct_id"]]
        struct_chain = row[defaults.DF_COLUMNS["struct_chain"]]
        uniprot_segments = row[defaults.DF_COLUMNS["uniprot_segments"]]
        struct_segments = row[defaults.DF_COLUMNS["struct_segments"]]

        #-------------------------------------------------------------#

        # If the structure is a homolog
        if source == "homologs":

            # The hit IDs will be the index of the data frame
            hit_id = index

        #-------------------------------------------------------------#

        # If data for the current structure are not available yet (it
        # can happen to encounter the same ID because of multiple
        # chains than can represent the protein of interest)
        if not struct_id in struct_data_dict.keys():

            # If the structure comes from the Protein Data Bank
            if source in ("structures", "homologs"):

                # Try to fetch the structure's data. Use only the first
                # model for multi-model structures
                try:
                    
                    struct_data = \
                        pdb_tools.get_pdb_data(struct_id = struct_id,
                                               model = 1)
                
                # If something went wrong
                except Exception as e:

                    # Log a warning
                    warnstr = \
                        "It was not possible to retrieve data for " \
                        f"the structure with ID '{struct_id}' from " \
                        "the Protein Data Bank. This structure " \
                        f"will be ignored. Error: {e}"
                    logger.warning(warnstr)
                    
                    # Go to the next structure
                    continue

                # Inform the user that data were successfully retrieved
                infostr = \
                    f"Data for the structure with ID '{struct_id}' " \
                    "were successfully retrieved from the Protein " \
                    "Data Bank."
                logger.info(infostr)

            #---------------------------------------------------------#

            # If the structure comes from the AlphaFold Protein
            # Structure Database
            elif source == "alphafold":

                # Try to fetch the structure's data
                try:
                    
                    struct_data = \
                        alphafold_tools.get_alphafold_data(\
                            struct_id = struct_id)
                
                # If something went wrong
                except Exception as e:

                    # Log a warning
                    warnstr = \
                        "It was not possible to retrieve data for " \
                        f"the structure with ID '{struct_id}' from " \
                        "the AlphaFold Protein Structure Database. " \
                        "This structure will be ignored. " \
                        f"Error: {e}"
                    logger.warning(warnstr)
                    
                    # Go to the next structure
                    continue

                # Inform the user that data were successfully retrieved
                infostr = \
                    f"Data for the structure with ID '{struct_id}' " \
                    "were successfully retrieved for the AlphaFold " \
                    "Protein Structure Database."
                logger.info(infostr)

            #---------------------------------------------------------#

            # Add data to the dictionary of structures' data
            struct_data_dict[struct_id] = struct_data

        #-------------------------------------------------------------#
        
        # For each SLiM still not found in a structure
        for slim in slims_to_find:

            #---------------------------------------------------------#

            # Get the SLiM's type
            s_type = slim[0]

            # Get where the SLiM starts
            s_start = slim[1]

            # Get the SLiM's sequence
            s_seq = slim[2]

            # Get where the SLiM ends
            s_end = slim[3]

            # If there is one more field
            if len(slim) == 5:
                
                # Get the SLiM's residues on which solvent
                # accessibility will be calculated
                s_acc_res = slim[4]
            
            # Otherwise
            else:
                
                # Set the residues to None
                s_acc_res = [None, None]

            #---------------------------------------------------------#
            
            # Try to find out if the SLiM is in the protein segments
            # covered by the structure
            try:
                
                slim_in_uniprot_segments = \
                    region_in_protein_segments(\
                        seq_range = (s_start, s_end),
                        uniprot_segments = uniprot_segments,
                        struct_segments = struct_segments)

            # If something went wrong
            except Exception as e:

                # Log a warning
                warnstr = \
                    "It was not possible to check if " \
                    f"{s_start}-{s_seq}-{s_end} is in the segments " \
                    "the UniProt sequence covered in the structure " \
                    f"with ID '{struct_id}'. This SLiM will be " \
                    f"ignored for this structure. Error: {e}"
                logger.warning(warnstr)
                
                # Go to the next SLiM
                continue

            #---------------------------------------------------------#

            # If the SLiM is in the structure
            if slim_in_uniprot_segments:

                # If the structure is of the protein of interest
                if source in ("structures", "alphafold"):
                    
                    # No need to check the structure for missing
                    # residues, since all segments with problems
                    # in the mapping were already excluded

                    slims_in_struct.append([s_type, \
                                           (s_start, s_seq, s_end),
                                            s_acc_res[0]])

                # If the structure is of a homolog
                elif source == "homologs":
                
                    # Try to find out if the sequence in the homolog
                    # corresponding to the SLiM in the protein
                    # of interest is in the structure
                    try:
                        
                        slim_in_struct_segments, start_in_t, \
                            seq_in_t, end_in_t = \
                                region_in_homolog_segments(\
                                    pdb_data = pdb_data,
                                    struct_id = struct_id,
                                    struct_chain = struct_chain,
                                    hit_al = als[hit_id],
                                    seq_range = (start, end))

                    # If something went wrong
                    except Exception as e:

                        # Log a warning
                        warnstr = \
                            "It was not possible to check whether " \
                            "the sequence " \
                            f"{s_start}-{s_seq}-{s_end} is in the " \
                            f"structure with ID '{struct_id}'. " \
                            "This SLiM will be ignored for this " \
                            f"structure. Error: {e}"
                        logger.warning(warnstr)
                        
                        # Go to the next SLiM
                        continue

                    # If the SLiM is in the structure
                    if slim_in_struct_segments:
                    
                        # Get the relative positions of the residues
                        # for which solvent accessibility will be
                        # calculated (with respect to the SLiM)
                        s_acc_res_rel = s_acc_res[1]

                        # If the positions are not None
                        if s_acc_res_rel is not None:

                            # Calculate the absolute positions of the
                            # residues for which solvent accessibility 
                            # will be calculated in the homolog
                            s_acc_res_abs = \
                                [start_in_t + res for res \
                                 in s_acc_res_rel]

                        # Otherwise
                        else:

                            # Set it to None
                            s_acc_res_abs = None

                        # Add the SLiM and the sequence corresponding
                        # to the SLiM in the homolog to the list of
                        # SLiMs found in the current structure
                        slims_in_struct.append(\
                            [s_type, \
                             (start_in_t, seq_in_t, end_in_t), \
                             (s_start, s_seq, s_end),
                             s_acc_res_abs])

        #-------------------------------------------------------------#

        # If no SLiMs are found in the structure
        if not slims_in_struct:

            # Inform the user of it
            infostr = \
                "No SLiMs found in the structure with ID " \
                f"'{struct_id}'."
            logger.info(infostr)
            
            # Go on to the next structure
            continue

        #-------------------------------------------------------------#

        # Inform the user of the SLiMs found in the current
        # structure
        slims_str = \
            ", ".join(\
                [f"{s_type}_{s_start}-{s_seq}-{s_end}" \
                 for (s_type, (s_start, s_seq, s_end), _) \
                 in slims_in_struct])
        infostr = \
            "The following SLiMs were found in the structure with " \
            f"ID '{struct_id}': {slims_str}."
        logger.info(infostr)

        #-------------------------------------------------------------#

        # Write out the file corresponding to the structure (it
        # is written out only if at least a SLiM has been found
        # in it)
        os.makedirs(struct_dir, exist_ok = True)
        struct_file = os.path.join(struct_dir, f"{struct_id}.pdb")
        with open(struct_file, "w") as out:
            out.write("\n".join(struct_data))

        # Inform the user that the structure was saved successfully
        infostr = \
            f"The structure with ID '{struct_id}' was successfully " \
            f"written in '{struct_file}'."
        logger.info(infostr)

        #-------------------------------------------------------------#

        # Add useful data to the list containing data about
        # the structure and the SLiMs found in it
        structures_and_slims.append((source, 
                                     struct_id,
                                     struct_chain,
                                     struct_file,
                                     slims_in_struct))

    #-----------------------------------------------------------------#

    # Return the list of PDB structures and the associated SLiMs
    return structures_and_slims, \
           [s for s in slims_to_find]



#######################################################################


def get_dssp_data(pdb_file,
                  dssp_params,
                  dssp_mapping,
                  wd):
    """Get data from DSSP's output file.

    Parameters
    ----------
    pdb_file : ``str``
        The PDB file of the structure on which DSSP should be run.

    dssp_params : ``dict``
        A dictionary containing the options to run DSSP.

    dssp_mapping : ``dict``
        A dictionary mapping symbols used by DSSP to identify
        single secondary structure elements to the same set
        of symbols, in case some of the secondary structure
        elements should be considered as belonging to the
        same class.

    wd : ``str``
        The working directory. DSSP's output will be written here.

    Returns
    -------
    df : ``pandas.DataFrame``
        A data frame containing the data parsed from DSSP's
        output file.

    fasta : ``dict``
        A dictionary mapping each chain found in DSSP's output
        file to its FASTA sequence.

    seqdssp : ``dict``
        A dictionary mapping each chain found in DSSP's output
        file to its FASTA-like sequence of DSSP symbols.
    """

    # Try to run DSSP
    try:
        
        _, dssp_output = \
            run_tools.run_dssp(pdb_file = pdb_file,
                               dssp_params = dssp_params,
                               wd = wd)

    # If something went wrong
    except Exception as e:

        # Set the error message
        errstr = \
            "It was not possible to run DSSP on the PDB file " \
            f"'{pdb_file}'."
        
        # Re-raise the exception
        raise Exception(errstr) from e

    # Inform the user that the DSSP run was successful
    infostr = \
        f"DSSP was successfully run on the PDB file '{pdb_file}'."
    logger.info(infostr)

    #-----------------------------------------------------------------#

    # Create the parser for the output of DSSP
    parser = dssp_tools.DSSPParser(dssp_mapping = dssp_mapping)

    # Try to parse the output
    try:
        
        dssp_data = parser.parse(dssp_output = dssp_output)

    # If something went wrong
    except Exception as e:

        # Set the error message
        errstr = \
            "It was not possible to parse data from DSSP's output " \
            f"file '{dssp_output}'."
        
        # Re-raise the exception
        raise Exception(errstr) from e

    # Inform the user that the output was successfully parsed
    infostr = \
        "Data were successfully parsed from DSSP's output file " \
        f"'{dssp_output}'." 
    logger.info(infostr)

    #-----------------------------------------------------------------#

    # Return the parsed data
    return dssp_data


def is_region_disordered(seq,
                         struct_id,
                         struct_chain,
                         fasta_dict,
                         seqdssp_dict,
                         min_len_disorder,
                         disorder_symbols):
    """Check if a portion of a protein sequence is in a disordered
    region, using data parsed from DSSP's output file.

    Parameters
    ----------
    seq : ``str``
        The portion of the protein sequence to be investigated.

    struct_id : ``str``
        The ID of the structure the DSSP data come from.

    struct_chain : ``str``
        The ID of the chain we are considering.

    fasta_dict : ``dict``
        A dictionary mapping each chain found in the structure
        on which DSSP was run to the chain's sequence.

        It is obtained when parsing DSSP's output file with
        ``dssp_tools.DSSPParser``.

    seqdssp_dict : ``dict``
        A dictionary mapping each chain found in the structure
        on which DSSP was run to the chain's FASTA-like sequence
        of DSSP symbols.

        It is obtained when parsing DSSP's output file with
        ``dssp_tools.DSSPParser``.

    min_len_disorder : ``int``
        The minimum length a disordered region should have to
        be considered a 'real' disordered region.

    disorder_symbols : ``list``
        A list of symbols representing disordered stated in
        the sequences provided in ``seqdssp_dict``.

    Returns
    -------
    is_disordered : ``bool`` or ``None``
        Whether the region is disordered.

        It is ``None`` if:
        - ``seq`` was not found in ``fasta_dict``.
        - Multiple occurrences of ``seq`` were found in ``fasta_dict``.
    """
    
    # Get the FASTA sequence of the chain of interest
    chain_f = fasta_dict[struct_chain]

    # Get the FASTA-like DSSP sequence of the chain of interest
    chain_d = seqdssp_dict[struct_chain]

    #-----------------------------------------------------------------#

    # Get all the occurrences of the sequence of interest in the FASTA
    # sequence of the chain
    occur = \
        [(i, chain_f[i:i+len(seq)], chain_f[i:i+len(seq)].find(seq)) \
         for i in range(len(chain_f)) \
         if chain_f[i:i+len(seq)].find(seq) != -1]

    #-----------------------------------------------------------------#

    # If no occurrence was found
    if len(occur) == 0:

        # Warn the user
        warnstr = \
            f"The sequence '{seq}' was not found in the structure " \
            f"with ID '{struct_id}', chain '{struct_chain}'."
        logger.warning(warnstr)

        # Return
        return

    #-----------------------------------------------------------------#
    
    # If multiple occurrences of the sequence were found
    elif len(occur) > 1:

        # Warn the user
        warnstr = \
            f"Multiple occurrences of '{seq}' were found in " \
            f"the structure with ID '{struct_id}', chain " \
            f"'{struct_chain}', and SLiMfast has no way of " \
            "determining which one is the one of interest."
        logger.warning(warnstr)

        # Return
        return

    #-----------------------------------------------------------------#

    # Get the starting position of the only occurrence found
    start_ix = occur[0][0]
    
    # Get the ending position of the only occurrence found
    end_ix = start_ix + len(seq)

    #-----------------------------------------------------------------#

    # Get the FASTA-like DSSP sequence of the occurrence
    seq_seqdssp = chain_d[start_ix:end_ix]

    #-----------------------------------------------------------------#

    # We want to be sure that the occurrence is found in a disordered
    # region at least 'min_len_disorder' residues long

    # Log the sequence we are working with now
    infostr = \
        f"Now working with sequence {start_ix+1}-{seq}-{end_ix} " \
        f"(DSSP sequence: '{seq_seqdssp}')..."
    logger.info(infostr)

    # Subtract the length of the sequence of interest from the
    # minimum length of the disordered region in order to obtain the
    # size of the extension of the sequence on both sides
    side_ext = math.ceil((min_len_disorder-len(seq)) / 2)
    
    # If the result is a negative extension (= the minimum length of
    # the disordered region is lower than the length of the sequence
    # of interest)
    if side_ext < 0:
        
        # The extension will be of length 0
        side_ext = 0

    #-----------------------------------------------------------------#

    # Get the putative starting position of the extended sequence of
    # interest
    ext_start_ix = start_ix - side_ext

    # Get the putative ending position of the extended sequence of
    # interest
    ext_end_ix = end_ix + side_ext

    #-----------------------------------------------------------------#

    # Try to extend the sequence on the left side
    try:
        
        _ext_seq_seqdssp = chain_d[ext_start_ix]
    
    # If there are not enough residues to extend the sequence on the
    # left side by the required number of residues, just extend it
    # as much as possible
    except IndexError:

        # Inform the user of the extension
        warnstr = \
            f"The sequence '{seq}' is close to the N-terminus of " \
            f"the structure with ID '{struct_id}', chain " \
            f"'{struct_chain}'. Therefore, the context on the left " \
            f"side is shorter than {side_ext} residues. The " \
            "available context will be considered."
        logger.warning(warnstr)
        
        # The new starting point will be the N-terminus
        ext_start_ix = 0

    #-----------------------------------------------------------------#

    # Try to extend the sequence on the right side
    try:
        
        _ext_seq_seqdssp = chain_d[ext_end_ix]
    
    # If there are not enough residues to extend the sequence on the
    # right side by the required number of residues, just extend it
    # as much as possible
    except IndexError:

        # Inform the user of the extension
        warnstr = \
            f"The sequence '{seq}' is close to the C-terminus of " \
            f"the structure with ID '{struct_id}', chain " \
            f"'{struct_chain}'. Therefore, the context on the " \
            f"right side is shorter than {side_ext} residues. " \
            "The available context will be considered."
        logger.warning(warnstr)
        
        # The new ending point will be the C-terminus
        ext_end_ix = len(chain_d) - 1

    #-----------------------------------------------------------------#

    # Get the extended FASTA sequence
    ext_seq = chain_f[ext_start_ix:ext_end_ix]

    # Get the extended DSSP sequence
    ext_seq_seqdssp = chain_d[ext_start_ix:ext_end_ix]
    
    # Log the extended sequence
    infostr = \
        f"Extended sequence (if 'min_len_disorder' > " \
        f"length_of_sequence): " \
        f"{ext_start_ix+1}-{ext_seq}-{ext_end_ix}. " \
        f"(DSSP sequence: '{ext_seq_seqdssp}')."
    logger.info(infostr)

    #-----------------------------------------------------------------#

    # If all residues in the extended sequence are disordered
    if set(ext_seq_seqdssp).issubset(set(disorder_symbols)):

        # Return True (the region is disordered)
        return True
    
    # Otherwise
    else:

        # Return False (the region is not disordered)
        return False


#######################################################################


def get_naccess_data(pdb_file,
                     naccess_params,
                     wd):
    """Get data from NACCESS' output file.

    Parameters
    ----------
    pdb_file : ``str``
        The PDB file of the structure on which NACCESS should be run.

    dssp_params : ``dict``
        A dictionary containing the options to run NACCCESS.

    wd : ``str``
        The working directory. NACCESS' output will be written here.

    Returns
    -------
    naccess_data : ``pandas.DataFrame``
        A data frame containing the data parsed from NACCESS'
        output file.
    """

    # Try to run NACCESS
    try:
        
        _, (naccess_output_abs, naccess_output_rel) = \
            run_tools.run_naccess(pdb_file = pdb_file,
                                  naccess_params = naccess_params,
                                  wd = wd)

    # If something went wrong
    except Exception as e:

        # Set the error message
        errstr = \
            "It was not possible to run NACCESS on the PDB file " \
            f"'{pdb_file}'."
        
        # Re-raise the exception
        raise Exception(errstr) from e

    #-----------------------------------------------------------------#

    # Create the parser for the output of NACCESS
    parser = naccess_tools.NACCESSParser()

    #-----------------------------------------------------------------#

    # Try to parse the output
    try:
        
        naccess_data = \
            parser.parse_rsa(naccess_output_rsa = naccess_output_rel)

    # If something went wrong
    except Exception as e:

        # Set the error message

        errstr = \
            "It was not possible to parse data from NACCESS' output " \
            f"file '{naccess_output_rel}'"
        
        # Re-raise the exception
        raise Exception(errstr) from e

    # Inform the user that the output was successfully parsed
    infostr = \
        "Data were successfully parsed from NACCESS' output file " \
        f"'{naccess_output_rel}'." 
    logger.info(infostr)

    #-----------------------------------------------------------------#

    # Return the parsed data
    return naccess_data


def get_naccess_region(model,
                       struct_chain,
                       naccess_data,
                       naccess_fields,
                       residues = None,
                       start = None,
                       end = None):
    """Get NACCESS data for a specific region and return them
    as a data frame.

    Parameters
    ----------
    model : ``int``
        The number of the model of interest in the structure.

    struct_chain : ``str``
        The ID of the chain where the region we are interested in lies.

    naccess_data : ``pandas.DataFrame``
        A data frame containing the data parsed from NACCESS' .rsa
        output file.

    naccess_fields : ``list``
        A list of strings representing the fields (= columns in
        the ``naccess_data`` data frame) that need to be reported
        in the data frame returned.

    residues : ``list``, optional
        A list of integers representing the sequence positions of
        the residues in the region of interest. 

    start : ``int``, optional
        The starting point of the sequence of interest.

        It is ignored if ``residues`` is provided.

    end : ``int``, optional
        The ending point of the sequence of interest.

        It is ignored if ``residues`` is provided.

    Returns
    -------
    naccess_region : ``pandas.DataFrame``
        A data frame containing the NACCESS data for the
        region of interest.
    """

    # If the user passed a list of residues of interest
    if residues is not None:

        # Select only those residues
        selection = \
            (naccess_data["model"] == model) & \
            (naccess_data["chain_id"] == struct_chain) & \
            (naccess_data["res_num"].isin(residues))

    #-----------------------------------------------------------------#

    # If the user passed both the starting and ending points of the
    # region of interest
    elif start is not None and end is not None:

        # Select only residues inside the region
        selection = \
            (naccess_data["model"] == model) & \
            (naccess_data["chain_id"] == struct_chain) & \
            (naccess_data["res_num"].isin(range(start, end+1)))

    #-----------------------------------------------------------------#

    # If the user passed only the starting position of the region of
    # interest
    elif start is not None and end is None:

        # Select all residues from the starting point onwards
        selection = \
            (naccess_data["model"] == model) & \
            (naccess_data["chain_id"] == struct_chain) & \
            (naccess_data["res_num"] > start)

    #-----------------------------------------------------------------#

    # If the user passed only the ending point of the region of
    # interest
    elif start is None and end is not None:

        # Select all residues before the ending point plus the ending
        # point itself
        selection = \
            (naccess_data["model"] == model) & \
            (naccess_data["chain_id"] == struct_chain) & \
            (naccess_data["res_num"] <= end)

    #-----------------------------------------------------------------#

    # Operate the selection
    naccess_region = naccess_data[selection][naccess_fields]

    #-----------------------------------------------------------------#

    naccess_region_dict = \
        {col : \
            "_".join([f"{index.split('_')[-1]}:{v}" \
                      for index, v in data.items()]) \
         for col, data in naccess_region.to_dict().items()}

    #-----------------------------------------------------------------#

    # Return the data for the region of interest
    return naccess_region_dict



#######################################################################


def get_results(structures_and_slims,
                wd,
                dssp_params,
                dssp_mapping,
                min_len_disorder,
                disorder_symbols,
                naccess_params = None):
    """Gather and return the results for the SLiMs found in
    specific structures.
    
    Parameters
    ----------
    structures_and_slims : ``list``
        A list of lists containing as many entries as the unique
        structures and associated unique chains in which SLiMs
        were found.

        Each entry consists of a list containing the following
        elements:

        - A string represented the type of 'source' represented
          by the structure. This can be ``"structures"``,
          ``"alphafold"``, or ``"homologs"`` depending on whether
          the structure is an experimental structure of the
          protein of interest from the Protein Data Bank, a
          structure of the protein of interest from the AlphaFold
          Protein Structure Database,  or the experimental structure
          of a homolog of the protein of interest from the Protein
          Data Bank.

        - A string representing the ID of the structure. This
          is the PDB ID for experimental structures coming from
          the Protein Data Bank and the AlphaFold ID for
          structures coming from the AlphaFold Protein Structure
          Database.

        - A string representing the chain ID of the structure's
          chain the SLiMs were found in.

        - A string representing the path to the PDB file containing
          the structure.

        - A list representing the SLiMs found in the given
          structure's chain.

    wd : ``str``
        The working directory. The output files for DSSP and,
        possibly, NACCESS, will be saved here.

    dssp_params : ``dict``
        A dictionary containing the options to run DSSP.

    dssp_mapping : ``dict``
        A dictionary mapping symbols used by DSSP to identify
        single secondary structure elements to the same set
        of symbols, in case some of the secondary structure
        elements should be considered as belonging to the
        same class.

    min_len_disorder : ``int``
        The minimum length (in number of residues) a disordered
        region must have to be considered disordered.

    disorder_symbols : ``list``
        A list of symbols belonging to the set used by DSSP to
        idenfity single secondary structure elements that should
        be considered symbols representing disordered regions.

    naccess_params :  ``dict``, optional
        A dictionary containing the options to run NACCESS.

        If not provided, NACCESS will not be run.

    Returns
    -------
    results_protein_df : ``pandas.DataFrame``
        A data frame containing the results for the SLiMs found
        in the experimental structures of the protein of interest
        coming from the Protein Data Bank.

    results_alphafold_df : ``pandas.DataFrame``
        A data frame containing the results for the SLiMs found
        in the structures of the protein of interest coming from
        the AlphaFold Protein Structure Database.

    results_homologs_df : ``pandas.DataFrame``
        A data frame containing the results for the SLiMs found
        in the structures of the homologs of the protein of interest
        coming from the Protein Data Bank.
    """

    # Create a dictionary to store the SLiMs found and annotated
    slims_results = \
        {"structures" : [], "alphafold" : [], "homologs" : []}

    #-----------------------------------------------------------------#

    # For each structure and SLiMs found in it
    for source, struct_id, struct_chain, pdb_file, slims_in_struct \
        in structures_and_slims:

        # Try to run DSSP, parse the DSSP output, and get the data
        try:
            _, fasta_dict, seqdssp_dict = \
                get_dssp_data(pdb_file = pdb_file,
                              dssp_params = dssp_params,
                              dssp_mapping = dssp_mapping,
                              wd = wd)

        # If something went wrong
        except Exception as e:

            # Log a warning
            warnstr = \
                "It was not possible to get the DSSP data for the " \
                f"structure with ID '{struct_id}'. This structure " \
                f"will be ignored. Exception: {e}"
            logger.warning(warnstr)
            
            # Go to the next structure
            continue

        #-------------------------------------------------------------#

        naccess_on_structs = set()

        #-------------------------------------------------------------#

        # For each SLiM found in the structure
        for slim_data in slims_in_struct:
    
            # If the SLiM comes from a structure of the protein
            # of interest
            if source in ("structures", "alphafold"):
                
                # Data will only be about the SLiM in the protein
                slim_name, (start, seq, end), slim_acc_res = slim_data
                start_q, seq_q, end_q = start, seq, end

            #---------------------------------------------------------#
                
            # If the SLiM comes from a structure of a homolog
            elif source == "homologs":
                
                # Data will be both about the SLiM in the protein
                # and about its corresponding sequence in the
                # homolog ("_q" stands for "query", therefore the
                # protein of interest)
                slim_name, (start, seq, end), \
                    (start_q, seq_q, end_q), slim_acc_res = slim_data

            #---------------------------------------------------------#

            # Try to check if the SLiM is in a disordered region
            try:
                
                is_disordered = \
                    is_region_disordered(\
                        seq = seq,
                        struct_id = struct_id,
                        struct_chain = struct_chain,
                        fasta_dict = fasta_dict,
                        seqdssp_dict = seqdssp_dict,
                        min_len_disorder = min_len_disorder,
                        disorder_symbols = disorder_symbols)

            # If something went wrong
            except Exception as e:

                # Log a warning
                warnstr = \
                    "It was not possible to check whether the " \
                    f"SLiM with UniProt sequence {start_q}-{seq_q}" \
                    f"-{end_q} is in a disordered region in the " \
                    f"structure with ID '{struct_id}' ({source}), " \
                    f"chain '{struct_chain}'. This SLiM will be " \
                    f"ignored. Exception: {e}"
                logger.warning(warnstr)
                
                # Go to the next structure
                continue

            #---------------------------------------------------------#

            # If the users passed residues on which NACCESS should be
            # run
            if slim_acc_res is not None:

                # Try to run NACCESS, parse the NACCESS output, and
                # get the data
                try:

                    naccess_data = \
                        get_naccess_data(pdb_file = pdb_file,
                                         naccess_params = naccess_params,
                                         wd = wd)

                # If something went wrong
                except Exception as e:

                    # Log a warning
                    warnstr = \
                        "It was not possible to get the NACCESS " \
                        "data for the structure with ID " \
                        f"'{struct_id}'. This structure will be " \
                        f"ignored. Exception: {e}"
                    logger.warning(warnstr)

                    # Go to the next structure
                    continue

                # Try to get NACCESS data for the region of interest
                try:

                    slim_acc = \
                        get_naccess_region(\
                            model = 1,
                            struct_chain = struct_chain,
                            residues = slim_acc_res,
                            naccess_data = naccess_data,
                            naccess_fields = defaults.NACCESS_COLS)

                # If something went wrong
                except Exception as e:

                    # Log a warning
                    warnstr = \
                        "It was not possible to get the NACCESS " \
                        "data for the SLiM with UniProt sequence " \
                        f"{start_q}-{seq_q}-{end_q} in the " \
                        f"structure with ID '{struct_id}' " \
                        f"('{source}'), chain '{struct_chain}'. " \
                        f"This SLiM will be ignored. Exception: {e}"
                    logger.warning(warnstr)

                    # Continue
                    continue

            # Otherwise
            else: 

                # The NACCESS data for the region will be an empty
                # dictionary
                slim_acc = {}
            
            #---------------------------------------------------------#

            # If the SLiM comes from a structure of the protein
            # of interest
            if source in ("structures", "alphafold"):
                
                # Save only information about the SLiM in the protein
                results = {"source" : source,
                           "struct_id" : struct_id,
                           "struct_chain" : struct_chain,
                           "slim_name" : slim_name,
                           "slim_start" : start_q,
                           "slim_sequence" : seq_q, 
                           "slim_end" : end_q,
                           "disordered_region" : str(is_disordered),
                           **slim_acc}

            #---------------------------------------------------------#

            # If the SLiM comes from a structure of a homolog
            elif source == "homologs":
                
                # Save information about both the SLiM in the
                # protein and its corresponding sequence in the
                # homolog
                results = {"source" : source,
                           "struct_id" : struct_id,
                           "struct_chain" : struct_chain,
                           "slim_name" : slim_name,
                           "slim_start" : start_q,
                           "slim_sequence" : seq_q, 
                           "slim_end" : end_q,
                           "slim_start_pdb" : start,
                           "slim_sequence_pdb" : seq,
                           "slim_end_pdb" : end,
                           "disordered_region" : str(is_disordered),
                           **slim_acc}

            #---------------------------------------------------------#

            # Add the SLiM to the results
            slims_results[source].append(results)

    #-----------------------------------------------------------------#

    # Create the data frames with the results
    results_protein_df = pd.DataFrame(slims_results["structures"])
    results_alphafold_df = pd.DataFrame(slims_results["alphafold"])
    results_homologs_df = pd.DataFrame(slims_results["homologs"])

    #-----------------------------------------------------------------#

    # Return the data frames
    return results_protein_df, \
           results_alphafold_df, \
           results_homologs_df


def get_disprot_annotations(uniprot_id,
                            slims_to_find):
    """Given a protein's UniProt ID, get the annotations provided
    by the DisProt database for a set of SLiMs of interest in the
    protein.

    Parameters
    ----------
    uniprot_id : ``str``
        The protein's UniProt ID.

    slims_to_find : ``list``
        A list containing the SLiMs for which DisProt annotations
        should be retrieved. Each SLiM should be represented by a
        list containing:

        - A string representing the SLiM's name.
        - An integer representing where the SLiM starts in the protein
          sequence.
        - A string representing the SLiM's sequence.
        - An integer representing where the SLiM ends in the protein
          sequence.
        - If a list of residues on which solvent accessibility should
          be calculated was passed, a tuple containing two tuples,
          the first containing a list of integers representing the
          absolute sequence positions of the residues, and the second
          containing a list of integers representing the relative
          sequence positions of the residues (with respect to the
          SLiM's sequence). If no list of residues was passed, the
          tuple will contain two ``None`` elements.

    Returns
    -------
    df_disprot_annot : ``pandas.DataFrame``
        A data frame containing the annotations retrieved. 
    """

    # Create the parser for DisProt annotations
    parser = disprot_tools.DisProtParser()

    #-----------------------------------------------------------------#

    # Try to get the annotations
    try:

        # Get the annotations for the given protein 
        df_disprot = parser.get_annotations(uniprot_id = uniprot_id)

        # Inform the user that the annotations were successfully
        # retrieved
        infostr = \
            "The DisProt annotations for the protein with UniProt " \
            f"ID '{uniprot_id}' were successfully retrieved."
        logger.info(infostr)

    # If something went wrong
    except Exception as e:

        # Warn the user
        errstr = \
            "It was not possbile to retrieve DisProt annotations " \
            f"for the protein with UniProt ID '{uniprot_id}'."
        logger.exception(errstr)

        # Re-raise the exception
        raise Exception(errstr) from e

    #-----------------------------------------------------------------#

    # If there were no annotations
    if df_disprot.empty:

        # Return the empty data frame
        return df_disprot

    #-----------------------------------------------------------------#

    # Get only the portion of the dataframe regarding disordered
    # regions, and only select the columns containing the starting
    # and ending points of the regions
    df_disprot_disorder = \
        df_disprot.loc[\
            df_disprot["term_name"] == "Disorder"][["start", "end"]]

    #-----------------------------------------------------------------#

    # Create an empty list to store the results
    results = []

    #-----------------------------------------------------------------#

    # For each SLiM
    for slim in slims_to_find:

        # Unpack the SLiM data
        s_type, s_start, s_seq, s_end, _ = slim
            
        # Try to find out if the SLiM is in one of the
        # disordered regions

        # Initialize the disordered flag to False
        is_disordered = False

        # For each disordered region
        for row in df_disprot_disorder.iterrows():
            
            # If the SLiM is included in the disordered region
            if start >= row["start"] and end <= row["end"]:

                # Turn the disordered flag on
                is_disordered = True

        # Save information about the SLiM
        results.append({"source" : "disprot",
                        "slim_type" : s_type,
                        "slim_start" : s_start,
                        "slim_sequence" : s_seq, 
                        "slim_end" : s_end,
                        "disordered_region" : str(is_disordered)})

    #-----------------------------------------------------------------#

    # Convert the results into a data frame and return it
    return pd.DataFrame(results)
