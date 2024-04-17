#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    disprot_tools.py
#
#    Utilities to retrieve, write and parse DisProt annotations.
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


#######################################################################


# Set the module's description
__doc__ = \
    "Utilities to retrieve, write and parse DisProt annotations."


#######################################################################


# Standard library
import json
import logging as log
import requests
# Third-party packages
import pandas as pd


#######################################################################


# Get the module's logger
logger = log.getLogger(__name__)


#######################################################################


class DisProtParser:

    """Class implementing a parser for DisProt annotations.
    """

    # The generic URL where the DisProt annotations are stored
    DISPROT_URL = "https://www.disprot.org/api/{:s}"

    # The names of the olumns of the data frame containing
    # the DisProt annotaions
    DF_COLS = \
        ["region_id", "start", "end",
         "term_id", "term_name", "term_namespace",
         "term_ontology",
         "ec_id", "ec_name", "ec_ontology",
         "reference_id", "reference_source",
         "curator_id", "curator_name", "curator_orcid",
         "date", "released", "cross_refs", "version"]


    def __init__(self):
        """Initialize the parser.
        """
        
        pass


    def get_annotations(self,
                        uniprot_id,
                        disprot_json_handle = None):
        """Get a data frame containing the DisProt annotations for
        a protein of interest, given the protein's UniProt ID.

        The annotations are directly parsed from the data retrieved
        from the DisProt database.

        Write a JSON file with the annotations in the original format,
        if requested.

        Parameters
        ----------
        uniprot_id : ``str``
            The protein's UniProt ID.

        disprot_json_handle : ``file_handle``
            The handle of the JSON file where to write the DisProt
            annotations in the original format. If not provided,
            the JSON file will not be written.

        Returns
        -------
        df_annot : ``pandas.DataFrame``
            The data frame with the annotations retrieved from
            DisProt.
        """

        # Create an empty dictionary to store the data about all
        # annotated regions of the protein (this dictionary will
        # be converted into a data frame)
        data = {}

        #-------------------------------------------------------------#

        # Get data from the server and convert the text into
        # a dictionary using JSON
        raw_data = \
            json.loads(requests.get(\
                self.DISPROT_URL.format(uniprot_id)).text)

        #-------------------------------------------------------------#

        # If there are data about the status of the entry
        if "status" in raw_data.keys():
            
            # If there are no annotations available yet
            if raw_data["status"] == "on-hold":
                
                # Warn the user
                warnstr = \
                    "The DisProt annotations for the protein with " \
                    f"UniProt ID '{uniprot_id}' have were not " \
                    "released yet (status: 'on-hold')."
                logger.warning(warnstr)
                
                # Return an empty data frame
                return pd.DataFrame()

            #---------------------------------------------------------#
            
            # If annotations do not exist
            elif raw_data["status"] == "not-exist":

                # Warn the user
                warnstr = \
                    "No DisProt annotations exist for the protein " \
                    f"with UniProt ID '{uniprot_id}' (status: " \
                    "'not-exist')."
                logger.warning(warnstr)
                
                # Return an empty data frame
                return pd.DataFrame()

        #-------------------------------------------------------------#
        
        # If the user requested the writing of a JSON file containing
        # the raw annotations
        if disprot_json_handle is not None:

            # Open the output file
            with disprot_json_handle as out:
                
                # Write out text data (JSON form)
                out.write(raw_data)

        #-------------------------------------------------------------#

        # Get data about the disordered regions
        regions = raw_data["regions"]

        #-------------------------------------------------------------#

        # For each region
        for num, region in enumerate(regions):

            # Create an empty dictionary to store data pertaining to
            # the current region
            region_data = {}

            # For each data field
            for field, value in region.items():

                # If the field contains either the HTML of the
                # reference or the statement
                if field == "reference_html" or field == "statement":

                    # Ignore it since it is a verbose field
                    continue

                # If the field holds the cross-references to other
                # databases
                elif field == "cross_refs":

                    # Format the cross-references
                    value = \
                        ":".join([f"{db['db']}={db['id']}" \
                                  for db in value])

                # Add the current field (and corresponding value)
                # to the dictionary storing data about the region
                region_data[field] = value

            # Add the dictionary to the dictionary containing data
            # for all regions
            data[num] = region_data

        #-------------------------------------------------------------#

        # Create the data frame and transpose it
        df = pd.DataFrame(data).T

        #-------------------------------------------------------------#

        # Get only those columns present in the data frame before
        # Q07817sorting them
        columns = \
            [col for col in self.DF_COLS if col in df.columns]

        #-------------------------------------------------------------#

        # Return the data frame with sorted columns
        return df[columns]
        