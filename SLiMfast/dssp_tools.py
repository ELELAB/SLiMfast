#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    dssp_tools.py
#
#    Contains a class implementing a parser to parse data from 
#    DSSP output files.
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
import logging as log
# Third-party packages
import pandas as pd



class DSSPParser:

    """Class implementing a parser for DSSP output files.
    """

    # Standard mapping of the DSSP symbols to the symbols
    # that will be used in the FASTA-like DSSP sequence
    # (default mapping is each symbol to itself)
    DSSP_MAPPING = \
        {"G" : "G", "I" : "I", "H" : "H", "B" : "B",
         "S" : "S", "E" : "E", "T" : "T"}

    # Rrepsentation of white spaces (= disordered residues)
    # in the FASTA-like DSSP sequence
    WHITESPACE_REPR = "-"


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

        # Get the mapping for DSSP symbols
        dssp_mapping = \
            dssp_mapping if dssp_mapping is not None \
            else self.DSSP_MAPPING

        # Get the whitespace representation
        whitespace_repr =  \
            whitespace_repr if whitespace_repr is not None \
            else self.WHITESPACE_REPR

        # Add the mapping for whitespaces to the final mapping that
        # will be used
        self.mapping = {**dssp_mapping, " " : whitespace_repr}


    def parse(self,
              dssp_output,
              verbose = False):
        """Parse the output of DSSP and return a dataframe
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
                    "Only parsing for non-verbose outputs has been " \
                    "implemented so far."
                raise NotImplementedError(errstr)

            # Create an empty list to store data from the DSSP file
            # that will be later converted into a dataframe
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

            for line in f:

                # Get rid of the newline delimiter
                line = line.rstrip("\n")

                # Start of the DSSP data
                if line[2] == "#":
                    # Turn on the boolean switch for parsing
                    start_parse = True
                    # Skip the current line
                    continue

                # If the boolean switch is off
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
                        f"{int(prev_res)+1} chain {prev_chain}."
                    log.warning(warnstr)
                    # Go to the next line
                    continue

                # Update the FASTA sequence
                tmp_fasta += aa

                # Update the FASTA-like DSSP sequence
                tmp_seqdssp += self.mapping[line[16]]

                # Append a dictionary of data found on this line
                # to the final list
                data.append(\
                    {"#" : line[1:5].lstrip(" "), 
                     "RESIDUE" : residue, 
                     "CHAIN" : chain,
                     "AA" : aa,
                     "STRUCTURE" : line[16:23].lstrip(" ").rstrip(" "),
                     "BP1" : line[24:28].lstrip(" "),
                     "BP2" : line[29:33].lstrip(" ").rstrip(" "),
                     "ACC" : line[34:37].lstrip(" "),
                     "NHO1" : line[38:49].lstrip(" "),
                     "OHN1" : line[50:61].lstrip(" "),
                     "NHO2" : line[62:73].lstrip(" "),
                     "OHN2" : line[74:85].lstrip(" "),
                     "TCO" : line[86:91].lstrip(" "),
                     "KAPPA" : line[92:97].lstrip(" "),
                     "ALPHA" : line[98:103].lstrip(" "),
                     "PHI" : line[104:109].lstrip(" "),
                     "PSI" : line[110:115].lstrip(" "),
                     "X-CA" : line[116:122].lstrip(" "),
                     "Y-CA" : line[123:129].lstrip(" "),
                     "Z-CA" : line[130:137].lstrip(" ")})

                # Update the holders of the previous chain
                # and previous residue
                prev_chain = chain
                prev_res = residue

            # Store the FASTA sequence and FASTA-like DSSP
            # sequence for the current chain
            fasta[chain] = tmp_fasta
            seqdssp[chain] = tmp_seqdssp

            # Return the dataframe and the two dictionaries
            return pd.DataFrame(data), fasta, seqdssp
