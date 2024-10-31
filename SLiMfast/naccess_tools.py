#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    naccess_tools.py
#
#    Utilities to parse data from NACCESS' output files.
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
__doc__ = "Utilities to parse data from NACCESS' output files."


#######################################################################


# Standard library
import logging as log
# Third-party packages
import pandas as pd
import re


#######################################################################


# Get the module's logger
logger = log.getLogger(__name__)


#######################################################################


class NACCESSParser:

    """Class implementing a parser for NACCESS' output files.

    The format of these files is described 
    `here <http://www.bioinf.manchester.ac.uk/naccess/>`_.
    """


    # A dictionary mapping residues' three-letter names to theis
    # one-letter names
    RESIDUES_3to1 = \
        {"ALA" : "A", "CYS" : "C", "ASP" : "D", "GLU" : "E",
         "PHE" : "F", "GLY" : "G", "HIS" : "H", "ILE" : "I",
         "LYS" : "K", "LEU" : "L", "MET" : "M", "ASN" : "N",
         "PRO" : "P", "GLN" : "Q", "ARG" : "R", "SER" : "S",
         "THR" : "T", "VAL" : "V", "TRP" : "W", "TYR" : "Y"}


    def parse_rsa(self,
                  naccess_output_rsa):
        """Parse the output files containing the per-residue
        absolute and relative solvent accessibility.

        Parameters
        ----------
        naccess_output_rsa : ``str``
            The .rsa output file produced by NACCESS.

        Returns
        -------
        df_rsa : ``pandas.DataFrame``
            A data frame containing the accessibilities parsed
            from the output file. 
        """

        # Open the file
        with open(naccess_output_rsa, "r") as f:

            # Set the number for the current model
            curr_model = 1

            #---------------------------------------------------------#

            # Initialize an empty dictionary to store the
            # data parsed from the file
            data = {}

            #---------------------------------------------------------#

            # For each line in the file
            for line in f:

                # Get the single items from the line
                items = [item for item in line.rstrip().split(" ") \
                         if item != ""]

                # If we are parsing a line containing residue-level 
                # information
                if items[0] == "RES":


                    # Handle the case where chain ID and residue number 
                    # are combined (e.g., 'A1200')
                    chain_res_combined = items[2]

                    
                    # Use a regex to separate the chain ID (one letter) 
                    # and the residue number (digits)
                    match = re.match(r"([A-Za-z])(\d+)",\
                                     chain_res_combined)

                    if match:

                        # This will be the chain ID, e.g., 'A'
                        chain_id = match.group(1)

                        # This will be the residue number, e.g., '1200'
                        res_num = match.group(2)

                        # Update the list to have the correct 
                        # list length
                        items[2] = chain_id
                        items.insert(3, res_num)

                    # In case of any unexpected format, 
                    # fallback to the original   
                    else:
                        chain_id = items[2]
                        res_num = items[3]

                    # Set the index for the current residue
                    index = f"{curr_model}_{chain_id}_\
                            {self.RESIDUES_3to1[items[1]]}_{res_num}"

                    # If it is already present in the data 
                    # (= we are parsing a new model)
                    if index in data:

                        # Update the number of the current model
                        curr_model += 1

                        # Re-define the index with the new model number
                        index = f"{curr_model}_{chain_id}_\
                                {self.RESIDUES_3to1[items[1]]}_\
                                {res_num}"

                    # Add the current residue to the dictionary
                    data[index] = {
                        "model": curr_model,
                        "chain_id": chain_id,
                        "res_num": int(res_num),
                        "res_name": items[1],
                        "abs_total_acc": float(items[4]),
                        "rel_total_acc": float(items[5]),
                        "abs_sc_total_acc": float(items[6]),
                        "rel_sc_total_acc": float(items[7]),
                        "abs_mc_total_acc": float(items[8]),
                        "rel_mc_total_acc": float(items[9]),
                        "abs_sc_nonpolar_acc": float(items[10]),
                        "rel_sc_nonpolar_acc": float(items[11]),
                        "abs_sc_polar_acc": float(items[12]),
                        "rel_sc_polar_acc": float(items[13]),
                    }
            #---------------------------------------------------------#

            # Convert the dictionary into a data frame and transpose it
            df_rsa = pd.DataFrame(data).T

            #---------------------------------------------------------#

            # Return the data frame
            return df_rsa