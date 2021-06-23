#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    disprot_tools.py
#
#    Contains a class implementing a parser to obtain, write and parse
#    DisProt annotations.
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
import requests
# Third-party packages
import pandas as pd



# Get the module logger
logger = log.getLogger(__name__)



class DisProtParser:

    """Class implementing a parser for DisProt annotations.
    """

    # Columns of the dataframe for the DisProt annotaions
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


    def get_annotations(uniprot_id,
                        disprot_json_handle = None):
        """Get a dataframe containing the DisProt annotations,
        parsed directly from the data retrieved from the DisProt
        database. Write a JSON file with the annotations in the
        original format, if requested.
        """

        # URL where DisProt data per UniProt ID are stored
        URL = "https://www.disprot.org/api/{:s}"

        # Get data from the server and convert the text into
        # a dictionary using JSON
        raw_data = \
            json.loads(requests.get(URL.format(uniprot_id)).text)
        
        # Whether to write out the annotations to a JSON
        # file
        if disprot_json_handle is not None:
            with disprot_json_handle as out:
                # Write out text data (JSON form)
                out.write(raw_data)

        # Get data about the disordered regions
        regions = raw_data["regions"]

        # Create an empty dictionary to store the data
        # about all regions (will be converted into a
        # dataframe)
        data = {}

        # For each region
        for num, region in enumerate(regions):

            # Create an empty dictionary to store data
            # pertaining to the current region
            region_data = {}

            # For each data field
            for field, value in region.items():

                # Ignore the HTML of the reference and the statement,
                # since they are verbose fields
                if field == "reference_html" or field == "statement":
                    continue

                # Format the cross-references to other databases
                elif field == "cross_refs":
                    value = \
                        ":".join([f"{db['db']}={db['id']}" \
                                  for db in value])

                # Add the current field (and corresponding value)
                # to the dictionary storing data about the region
                region_data[field] = value

            # Add the dictionary to the dictionary containing
            # data for all regions
            data[num] = region_data

        # Create and return the dataframe
        return pd.DataFrame(data).T[self.DF_COLS]
