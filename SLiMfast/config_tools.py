#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    config_tools.py
#
#    Utilities to check configurations.
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


# Set the module's description
__doc__ = "Utilities to check configurations."


#######################################################################


# Standard library
import logging as log
# SLiMfast
from . import defaults


#######################################################################


# Get the module's logger
logger = log.getLogger(__name__)


########################## PRIVATE UTILITIES ##########################


def _merge_templates(templates):
    """Merge two or more templates into one.

    Parameters
    ----------
    templates : ``list``
        A list of configuration templates.

    Returns
    -------
    merged_template : ``dict``
        A template obtained by merging all the input templates.
    """

    # Initialize the merged template to the first template
    merged = templates[0]

    #-----------------------------------------------------------------#

    # Define the recursion
    def recurse(template_1,
                template_2):

        # Initialize the result to the first template provided
        result = dict(template_1)

        # For each key, value pair in the second template
        for key, value in template_2.items():

            # If the key is also present in the first template and
            # is associated with a dictionary in both templates
            if key in template_1 \
            and isinstance(template_1[key], dict) \
            and isinstance(value, dict):

                # Recursively merge the value
                result[key] = recurse(template_1 = template_1[key],
                                      template_2 = value)

            # Otherwise
            else:

                # Overwrite an existing key or add a key from the
                # second template
                result[key] = value

        # Return the resulting template
        return result

    #-----------------------------------------------------------------#

    # For each of the other templates
    for current_template in templates[1:]:

        # Iteratively create the merged dictionari
        merged = recurse(template_1 = merged,
                         template_2 = current_template)

    #-----------------------------------------------------------------#

    # Return the merged template
    return merged


########################## PUBLIC UTILITIES ###########################


def check_against_template(config,
                           config_type):
    """Check the configuration against the configuration's own
    template.

    Parameters
    ----------
    config : ``dict``
        The configuration loaded from a YAML configuration file.

    config_type : ``str``, {``"slimfast_get_structures"``,
        ``"slimfast_get_alphafold"``, ``"slimfast_get_homologs"``,
        ``"slimfast_find_slims"``, ``"slimfast_run_pipeline"``}
        The type of configuration we are checking.

    Returns
    -------
    config : ``dict``
        The configuration loaded from the file provided by
        the user, if all checks were successful, possibly updated.
    """

    # The options of type 'dict'
    DICT_OPTIONS = \
        [("structures", "alphafold", "df", "filter"),
         ("structures", "homolog", "df", "filter"),
         ("hhblits", "run", "options"),
         ("hhsearch", "run", "options"),
         ("dssp", "run", "options"),
         ("naccess", "run", "options")]

    # Non-mandatory options
    OPTIONAL_OPTIONS = \
        [("structures", "protein", "df", "sort"),
         ("structures", "alphafold", "df", "sort"),
         ("structures", "homolog", "df", "sort"),
         ("structures", "alphafold", "df", "filter"),
         ("structures", "homolog", "df", "filter")]

    #-----------------------------------------------------------------#

    # Define the recursion
    def recurse(config,
                config_copy,
                template,
                path = ()):

        # For each key, value pair in the template
        for key, value in template.items():

            # If the value is a dictionary
            if isinstance(value, dict) \
            and path + (key,) not in DICT_OPTIONS:

                # Set the string representing the path to the
                # section
                sectstr = ":".join(path + (key,))

                # If they key is not in the configuration
                if key not in config \
                and path + (key,) not in OPTIONAL_OPTIONS:

                    # Raise an error
                    errstr = \
                        f"The section '{sectstr}' must be present " \
                        "in the configuration."
                    raise KeyError(errstr)

                # Recurse through the configuration
                recurse(config = config[key],
                        config_copy  = config_copy[key],
                        template = value,
                        path = path + (key,))

            #---------------------------------------------------------#

            # Otherwise
            else:

                # Set the string representing the path to the
                # option
                optstr = ":".join(path + (key,))

                # If the value is a data type
                if isinstance(value, type):

                    # If the key is not in the configuration
                    if key not in config \
                    and path + (key,) not in OPTIONAL_OPTIONS:

                        # Raise an error
                        errstr = \
                            f"The option '{optstr}' must be " \
                            "present in the configuration."
                        raise KeyError(errstr)

                    # If the value associated with the key in the
                    # configuration is not of the right type
                    if not isinstance(config[key], value):
                        
                        # Raise an error
                        errstr = \
                            f"The option '{optstr}' must be of " \
                            f"{value}."
                        raise ValueError(errstr)

                # Otherwise
                else:

                    # Get the option's type and its default value
                    opt_type, opt_default = value

                    # If the key is not in the configuration
                    if key not in config:

                        # Update the copy of the configuration with
                        # the default value for the option
                        config_copy[key] = value

                    # Otherwise
                    else:

                        # If the value associated with the key in the
                        # configuration is not of the right type
                        if not isinstance(config[key], opt_type):
                            
                            # Raise an error
                            errstr = \
                                f"The option '{optstr}' must be of " \
                                f"{opt_type}."
                            raise ValueError(errstr)

        #-------------------------------------------------------------#

        # Return the copy of the configuration, possibly updated
        return config_copy

    #-----------------------------------------------------------------#

    # If the configuration is for 'slimfast_run_pipeline'
    if config_type == "slimfast_run_pipeline":

        # Create the template by merging the ones for the other
        # configurations
        template = \
            _merge_templates(\
                templates = list(defaults.CONFIG_TEMPLATES.values()))

    # Otherwise
    else:

        # Pick the relevant template
        template = defaults.CONFIG_TEMPLATES[config_type]

    #-----------------------------------------------------------------#

    # Recurse and return the result
    return recurse(config = config,
                   config_copy = config.copy(),
                   template = template)
