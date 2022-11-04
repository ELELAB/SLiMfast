#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    dssp_tools.py
#
#    Contains a class implementing a parser to parse data from 
#    DSSP output files.
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
import logging as log
# Third-party packages
import pandas as pd



# Get the module logger
logger = log.getLogger(__name__)



class DSSPParser:

    """Class implementing a parser for DSSP output files.
    """

    # Standard mapping of the DSSP symbols to the symbols
    # that will be used in the FASTA-like DSSP sequence
    # (the default mapping is each symbol to itself)
    DSSP_MAPPING = \
        {"G" : "G", "I" : "I", "H" : "H", "B" : "B",
         "S" : "S", "E" : "E", "T" : "T"}

    # Representation of white spaces (= disordered residues)
    # in the FASTA-like DSSP sequence
    WHITESPACE_REPR = "-"

    # Data types for each column of the final data frame
    DTYPES = {\
        "RESIDUE" : "int64",
        "CHAIN" : "object",
        "AA" : "object",
        "STRUCTURE" : "object",
        "3-TURNS_HELIX" : "object",
        "4-TURNS_HELIX" : "object",
        "5-TURNS_HELIX" : "object",
        "GEOMETRICAL_BEND" : "object",
        "CHIRALITY" : "object",
        "BETA_BRIDGE_LABEL" : "object",
        "BP1" : "int64",
        "BP2" : "object",
        "ACC" : "int64",
        "N-H-->O_1" : "object",
        "O-->H-N_1" : "object",
        "N-H-->O_2" : "object",
        "O-->H-N_2" : "object",
        "TCO" : "float64",
        "KAPPA" : "float64",
        "ALPHA" : "float64",
        "PHI" : "float64",
        "PSI" : "float64",
        "X-CA" : "float64",
        "Y-CA" : "float64",
        "Z-CA" : "float64"}


    def __init__(self,
                 dssp_mapping = None,
                 whitespace_repr = None):
        """Initialize the DSSP output parser.
        
        Parameters:
        -----------
        dssp_mapping : `dict` or `None`, default: `None`
            Dictionary mapping the 7 symbols of the DSSP dictionary
            to custom symbols.

        whitespace_repr : `str`, default: `"-"`
            How the whitespaces (assigned by DSSP to unstructured
            residues) are represented in the FASTA-like representation
            of the DSSP sequence.
        """

        # Get the mapping for DSSP symbols (if none provided, fall
        # back on the default)
        dssp_mapping = \
            dssp_mapping if dssp_mapping is not None \
            else self.DSSP_MAPPING

        # Get the whitespace representation (if none provided, fall
        # back on the default)
        whitespace_repr =  \
            whitespace_repr if whitespace_repr is not None \
            else self.WHITESPACE_REPR

        # Add the mapping for whitespaces to the final mapping that
        # will be used
        self.mapping = {**dssp_mapping, " " : whitespace_repr}


    def parse(self,
              dssp_output,
              verbose = False):
        """Parse the output of DSSP and return a data frame
        with the data.

        Parameters
        ----------
        dssp_output : `str`
            File where the DSSP output is stored.

        verbose : `bool`, default: `False`
            Whether the DSSP output has been generated in verbose
            mode.
        """

        with open(dssp_output, "r") as f:       

            # If 'verbose' has been set to True
            if verbose:
                
                # Raise an exception, since parsing for output files
                # generated in verbose mode has not been implemented
                # yet
                errstr = \
                    "Only parsing for non-verbose DSSP outputs has " \
                    "been implemented so far."
                raise NotImplementedError(errstr)

            # Create an empty list to store data from the DSSP file
            # that will be later converted into a data frame
            data = []

            # Create an empty dictionary to store the FASTA sequence
            # (per chain) of the structure on which DSSP has been run
            seqdssp = {}

            # Create an empty dictionary to store the FASTA-like DSSP
            # sequence (per chain) of the structure on which DSSP
            # has been run
            fasta = {}

            # Initialize the boolean switch that indicates when
            # to start parsing to False
            start_parse = False

            # Initialize the string that will temporarily hold the
            # FASTA sequence for a single chain to an empty string
            tmp_fasta = ""

            # Initialize the string that will temporarily hold the
            # FASTA-like DSSP sequence for a single chain to an
            # empty string 
            tmp_seqdssp = ""

            # Initialize the string that will temporarily hold the
            # ID of the chain to which the residue prior to the one
            # currently parsed belong to an empty string
            prev_chain = ""

            # Initialize the string that will temporarily hold the
            # residue prior to the one currently parsed to an
            # empty string
            prev_res = ""

            # For each line in the DSSP output
            for line in f:

                # Get rid of the newline delimiter
                line = line.rstrip("\n")

                # Start of the DSSP data
                if line[2] == "#":
                    
                    # Turn on the boolean switch for parsing
                    start_parse = True
                    
                    # Skip the current line
                    continue

                # If the boolean switch for starting parsing is off
                if not start_parse:
                    
                    # Skip the current line
                    continue

                # Get residue number, chain ID and residue type
                residue = line[6:10].lstrip(" ")
                chain = line[11].lstrip(" ")
                aa = line[13:15].rstrip(" ")

                # If you reached a chain break
                if aa == "!*":
                    
                    # Save the FASTA-like DSSP sequence of
                    # the chain that just ended
                    seqdssp[prev_chain] = tmp_seqdssp
                    
                    # Save the FASTA sequence of the chain
                    # that just ended
                    fasta[prev_chain] = tmp_fasta
                    
                    # Re-initialize the temporary holders of
                    # the FASTA sequence and of the FASTA-like
                    # DSSP sequence
                    tmp_fasta = ""
                    tmp_seqdssp = ""
                    
                    # Go to the next line
                    continue
                
                # If you found a missing residue
                if aa == "!":
                    
                    # Warn the user
                    warnstr = \
                        f"Missing residue at position " \
                        f"{int(prev_res)+1} chain {prev_chain} " \
                        f"in {dssp_output}."
                    logger.warning(warnstr)
                    
                    # Go to the next line
                    continue

                # Update the FASTA sequence
                tmp_fasta += aa

                # Update the FASTA-like DSSP sequence
                tmp_seqdssp += self.mapping[line[16]]

                # Append a dictionary of data found on this line
                # to the final list
                data.append(\
                    {"#" : int(line[1:5].lstrip(" ")), 
                     "RESIDUE" : int(residue), 
                     "CHAIN" : chain,
                     "AA" : aa,
                     "STRUCTURE" : line[16],
                     "3-TURNS_HELIX" : line[18],
                     "4-TURNS_HELIX" : line[19],
                     "5-TURNS_HELIX" : line[20],
                     "GEOMETRICAL_BEND" : line[21],
                     "CHIRALITY" : line[22],
                     "BETA_BRIDGE_LABEL" : line[23:25],
                     "BP1" : int(line[25:29].lstrip(" ")),
                     "BP2" : line[29:34].lstrip(" ").rstrip(" "),
                     "ACC" : int(line[34:38].lstrip(" ")),
                     "N-H-->O_1" : line[38:50].lstrip(" ").rstrip(" "),
                     "O-->H-N_1" : line[50:61].lstrip(" ").rstrip(" "),
                     "N-H-->O_2" : line[62:73].lstrip(" ").rstrip(" "),
                     "O-->H-N_2" : line[74:85].lstrip(" ").rstrip(" "),
                     "TCO" : float(line[85:91].lstrip(" ")),
                     "KAPPA" : float(line[91:97].lstrip(" ")),
                     "ALPHA" : float(line[97:103].lstrip(" ")),
                     "PHI" : float(line[103:109].lstrip(" ")),
                     "PSI" : float(line[109:115].lstrip(" ")),
                     "X-CA" : float(line[115:122].lstrip(" ")),
                     "Y-CA" : float(line[122:129].lstrip(" ")),
                     "Z-CA" : float(line[129:137].lstrip(" "))})

                # Update the holders of the previous chain
                # and previous residue
                prev_chain = chain
                prev_res = residue

            # Store the FASTA sequence and FASTA-like DSSP
            # sequence for the current chain
            fasta[chain] = tmp_fasta
            seqdssp[chain] = tmp_seqdssp

            # Create the data frame
            df = pd.DataFrame(data).set_index("#").astype(self.DTYPES)

            # Return the data frame and the two dictionaries
            return df, fasta, seqdssp
