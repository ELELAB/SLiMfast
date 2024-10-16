#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    hhseach_tools.py
#
#    Utilities to parse data from the output files produced by
#    HHsearch and manipulate them.
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
__doc__ = \
    "Utilities to parse data from the output files produced by " \
    "HHsearch manipulate them."


#######################################################################


# Standard library
import logging as log
import re
# Third-party packages
import pandas as pd
# SLiMfast
from . import defaults


#######################################################################


# Get the module logger
logger = log.getLogger(__name__)


#######################################################################


class HHsearchParser:

    """Class implementing a parser for the output file produced by
    HHsearch.
    """

    # The format for data provided in scientific notation
    SCI_FMT = "{:,.1E}"

    # The names of the columns in the hits' table
    COLNAMES_HITS = \
        ["num", "hit", "description", "prob", "e_value",
         "p_value", "score", "ss", "cols",
         defaults.DF_COLUMNS["uniprot_segments"], 
         defaults.DF_COLUMNS["struct_segments"]]

    # The names of the columns for the data retrieved from the
    # aligments
    COLNAMES_ALS = \
        ["identities", "similarity", "sum_probs", "template_neff",
         "organism"]

    # Column names for the additional columns that will contain the
    # structure's PDB ID and chain ID
    COLNAMES_STRUCT = \
        [defaults.DF_COLUMNS["struct_id"],
         defaults.DF_COLUMNS["struct_chain"]]

    # The name of the additional column that will contain information
    # about the sequence coverage
    COLNAME_COV = "uniprot_coverage_percent"

    # The separator to be used between the hit number and the hit ID
    # (hit number + hit ID will be the hit's unique identifier)
    HIT_NUM_ID_SEP = "_"

    # The data types for each column of the final data frame
    DTYPES = \
        {"num" : "int64",
         "hit" : "object",
         "description" : "object",
         "prob" : "float64",
         "e_value" : "object",
         "p_value" : "object",
         "score" : "float64",
         "ss" : "float64",
         "cols" : "int64",
         defaults.DF_COLUMNS["uniprot_segments"] : "object",
         defaults.DF_COLUMNS["struct_segments"] : "object",
         defaults.DF_COLUMNS["struct_id"] : "object",
         defaults.DF_COLUMNS["struct_chain"] : "object",
         "identities" : "float64",
         "similarity" : "float64",
         "sum_probs" : "float64",
         "template_neff" : "float64",
         "organism" : "object",
         COLNAME_COV : "float64"}


    def __init__(self):
        """Initialize the parser.
        """
        
        pass


    def _get_section_indexes(self,
                             f):
        """Get the indexes of the hits' and alignments' sections.

        Parameters
        ----------
        f : ``file handle``
            The file handle.

        Returns
        -------
        hits_indexes : ``tuple``
            The starting and ending point of the hits' section
            (line numbers).

        als_indexes : ``tuple``
            The starting and ending point of the alignments'
            section.
        """

        # Set a boolean switch to indicate when to parse the table
        # containing the hits; it is off at the beginning since
        # the table comes after the metadata
        parse_hits = False

        #-------------------------------------------------------------#
            
        # Set a boolean switch to indicate when to parse the
        # alignments; it is off at the beginning since the alignments
        # come after the hits' table
        parse_als = False

        #-------------------------------------------------------------#

        # Initialize the index of the line where the hits' table starts
        # to 0
        start_hits_ix = 0
        
        # Initialize the index of the line where the hits' table ends
        # to 0
        end_hits_ix = 0

        #-------------------------------------------------------------#
        
        # For each line
        for ix, line in enumerate(f):

            # If we encounter an empty line
            if re.match(r"^\s*$", line):

                # If it is the first empty line encountered (= we
                # have finished parsing the metadata)
                if not parse_hits and not parse_als:
                    
                    # Turn on the parsing for the hits' table since
                    # it is the next section
                    parse_hits = True
                    
                    # The starting index for the hits' table section
                    # is the one of the next line
                    start_hits_ix = ix + 1
                    
                    # Ignore the current (empty) line
                    continue

            #---------------------------------------------------------#

            # If we are parsing the hits' table section
            if parse_hits:

                # If we encounter a blank line
                if re.match(r"^\s*$", line):
                    
                    # Turn off the parsing of the hits' table
                    parse_hits = False
                    
                    # Turn on the parsing for the alignments
                    parse_als = True
                    
                    # The ending index for the hits' table is the one
                    # of the current line (the interval is open, so
                    # the table ends on the previous line)
                    end_hits_ix = ix
                    
                    # Ignore the current (empty) line
                    continue

        #-------------------------------------------------------------#

        # Return the indexes of the hits' and the aligments' sections
        return (start_hits_ix, end_hits_ix), (end_hits_ix+1, ix)


    def _parse_hits(self,
                    f,
                    range_indexes):
        """Parse the hits' section.

        Parameters
        ----------
        f : ``file handle``
            The file handle.

        range_indexes : ``tuple``
            A tuple containing the starting and ending point of the
            hits' section (line numbers).

        Returns
        -------
        df_hits : ``pandas.DataFrame``
            A data frame containing the data parsed from the hits'
            section.
        """

        # Set the name of the field that will contain the segments
        # covered by the structure using the UniProt numbering
        uniprot_segments_col = defaults.DF_COLUMNS["uniprot_segments"]

        # Set the name of the field that will contain the segments
        # covered by the structure using the structure's own
        # numbering
        struct_segments_col = defaults.DF_COLUMNS["struct_segments"]

        #-------------------------------------------------------------#

        # Create an empty dictionary to store data
        data = []

        #-------------------------------------------------------------#

        # Set a boolean switch for when to start parsing the hits'
        # section 
        parse = False

        #-------------------------------------------------------------#

        # Unpack the range in starting and ending index
        start_ix, end_ix = range_indexes

        #-------------------------------------------------------------#

        # Go to the beginning of the file
        f.seek(0)

        #-------------------------------------------------------------#
        
        # For each line
        for i, line in enumerate(f):

            # If we reached the beginning of the hits' section
            if i == start_ix + 1:

                # Turn on the switch for parsing
                parse = True

            #---------------------------------------------------------#

            # If we reached the end of the hits' section
            elif i == end_ix:

                # Exit the loop
                break

            #---------------------------------------------------------#

            # If the parsing switch is not active
            if not parse:

                # Continue
                continue

            #---------------------------------------------------------#

            # Split the line with whitespaces and take all non-empty
            # elements
            line = [i for i in line.split(" ") if i != ""]

            #---------------------------------------------------------#

            # Get the hit's ID
            hit = \
                line[1].split("_")[0].lower() \
                + "." \
                + line[1].split("_")[1]

            #---------------------------------------------------------#

            # If the structure's segments are stored in the last
            # elements of the line
            if "-" in line[-1]:

                # Add the data to the final list
                data.append(\
                    {"num" : int(line[0]),
                     "hit" : hit,
                     "description" : " ".join(line[2:-8]),
                     "prob" : float(line[-8]),
                     "e_value" : self.SCI_FMT.format(float(line[-7])),
                     "p_value" : self.SCI_FMT.format(float(line[-6])),
                     "score" : float(line[-5]),
                     "ss" : float(line[-4]),
                     "cols" : int(line[-3]),
                     uniprot_segments_col : line[-2],
                     struct_segments_col : line[-1].split("(")[0]})

            # If the structure's segments are stored in the
            # second-to-last element of the line
            else:

                # Add the data to the final list
                data.append(\
                    {"num" : int(line[0]),
                     "hit" : hit,
                     "description" : " ".join(line[2:-9]),
                     "prob" : float(line[-9]),
                     "e_value" : self.SCI_FMT.format(float(line[-8])),
                     "p_value" : self.SCI_FMT.format(float(line[-7])),
                     "score" : float(line[-6]),
                     "ss" : float(line[-5]),
                     "cols" : int(line[-4]),
                     uniprot_segments_col : line[-3],
                     struct_segments_col : line[-2]})

        #-------------------------------------------------------------#

        # Convert the dictionary containing the hits into a data frame
        df_hits = pd.DataFrame(data)

        #-------------------------------------------------------------#

        # Set the names of the columns
        df_hits.columns = self.COLNAMES_HITS

        #-------------------------------------------------------------#

        # Add two columns storing just the PDB ID and the chain
        df_hits[self.COLNAMES_STRUCT] = \
            df_hits["hit"].str.split(pat = ".",
                                     n = 1,
                                     expand = True)

        #-------------------------------------------------------------#

        # The rows' names are composed by both the hit number and the
        # structure's PDB ID (with chain) to both ensure their
        # uniqueness and be informative at the same time
        df_hits.index = \
            [str(n) + self.HIT_NUM_ID_SEP + h for n, h \
             in zip(df_hits["num"], df_hits["hit"])]

        #-------------------------------------------------------------#

        # Return the data frame
        return df_hits


    def _parse_als(self,
                   f,
                   range_indexes):
        """Parse the alignments' sections and return a dictionary of
        dictionaries containing the strings corresponding to
        each alignment. The alignment-to-dictionary transformation
        is shown below.

        Parameters
        ----------
        f : ``file_handle``
            The file handle.

        range_indexes : ``tuple``
            A tuple containing the starting and ending point of the
            alignments' section (line numbers).

        Returns
        -------
        alignments : ``dict``
            A dictionary of dictionaries containing the alignments.
        
        Notes
        -----
        How an alignment is the alignments' section gets converted
        into a dictionary.

        *Alignment*

        No 1
        > <target data>
        <other data>

        Q ss_pred             cCCCCCCCCCCceeeEEEEcCCCCceEEEEEC
        Q query             1 ATQSKPTPVKPNYALKFTLAGHTKAVSSVKFS  32 (32)
        Q Consensus         1 ~~~~~~~~~~~~~~~~~~~~~~~~~i~~~~~s  32 (32)
                              +++++.+....+...+.+.+|...|.+++|+|
        T Consensus         1 ~~~~~~~~~~~~~~~~~~~~~~~~~i~~~~~~  32 (32)
        T target            1 ATQSKPTPVKPNYALKFTLAGHTKAVSSVKFS  32 (32)
        T ss_dssp             --------CCCCCEEEEEECCCSSCEEEEEEC
        T ss_pred             CCCCCCCCCCCcEEEEeecCCcCCCeeEEEEc

        
        *Dictionary*
        
        {hit_id : 
            {"Q" : {"ss_pred" : "cCCCCCCCCCCceeeEEEEcCCCCceEEEEEC",
                    "query" : "ATQSKPTPVKPNYALKFTLAGHTKAVSSVKFS",
                    "Consensus" : "~~~~~~~~~~~~~~~~~~~~~~~~~i~~~~~s"},

             "quality" : "+++++.+....+...+.+.+|...|.+++|+|",
             
             'T' : {"Consensus" : "~~~~~~~~~~~~~~~~~~~~~~~~~i~~~~~~",
                    "template" : "ATQSKPTPVKPNYALKFTLAGHTKAVSSVKFS",
                    "ss_dssp" : "--------CCCCCEEEEEECCCSSCEEEEEEC",
                    "ss_pred" : "CCCCCCCCCCCcEEEEeecCCcCCCeeEEEEc"}
            }
        }
        """

        # Set a boolean switch indicating whether to parse the
        # alignments to False
        parse = False

        #-------------------------------------------------------------#

        # Set the starting and ending points for query and template
        # to None
        start_q = None
        start_t = None
        end_q = None
        end_t = None

        #-------------------------------------------------------------#

        # Set the hit ID to None
        hit_id = None

        #-------------------------------------------------------------#

        # Set the number of residues in the query to None
        num_res = None

        #-------------------------------------------------------------#

        # Initialize an empty dictionary to store the information on
        # the hits gathered from the alignments
        info_from_als = {}

        #-------------------------------------------------------------#

        # Create an empty dictionary to store the alignments
        als = {}

        #-------------------------------------------------------------#        

        # Go to the beginning of the file
        f.seek(0)

        #-------------------------------------------------------------#

        # For each line
        for ix, line in enumerate(f):

            # If we reached the beginning of the alignments' section,
            if ix == range_indexes[0]:

                # Turn on the boolean switch
                parse = True

            #---------------------------------------------------------#

            # If parsing is not active
            if not parse:

                # Continue (= skip all sections of the file that do
                # not contain the alignments)
                continue

            #---------------------------------------------------------#

            # If we encounter an empty line
            elif re.match(r"^\s*$", line) \
            and (not ix == range_indexes[1]):

                # Ignore it
                continue

            #---------------------------------------------------------#

            # If we just entered a new aligment or we reached the end
            # of the file
            elif line.startswith("No") \
            or (ix == range_indexes[1] and re.match(r"^\s*$", line)):

                # If there are data from a previous hit to be stored
                if hit_id is not None:

                    # If no info for the coverage of the UniProt
                    # sequence by the alignment was stored for the hit
                    if not self.COLNAME_COV in \
                        info_from_als[hit_id].keys() \
                    and num_res is not None:
                                    
                        # Get the coverage of the UniProt sequence by
                        # the current alignment in percent points
                        cov_percent = \
                            float("{:.1f}".format(\
                                (int(aligned_cols) / num_res) * 100))

                        # Save this info for the current alignment
                        info_from_als[hit_id].update(\
                            {self.COLNAME_COV : cov_percent})

                    #-------------------------------------------------#
                    
                    # Save the start and end numbers of both the query
                    # and the template/target
                    als[hit_id]["Q"]["start"] = start_q
                    als[hit_id]["Q"]["end"] = end_q
                    als[hit_id]["T"]["start"] = start_t
                    als[hit_id]["T"]["end"] = end_t

                    #-------------------------------------------------#

                    # If we have reached the end of the file
                    if ix == range_indexes[1] \
                    and re.match(r"^\s*$", line):

                        # Break out of the loop (there is no need to
                        # run the next lines of code in this "if"
                        # block since, if we up here, it is because
                        # we only have one entry in our results)
                        break

                    #-------------------------------------------------#

                    # Re-initialize all the starting and ending points
                    # to None for the next hit
                    start_q = None
                    end_q = None
                    start_t = None
                    end_t = None

                #-----------------------------------------------------#
                
                # Get the hit number
                hit_num = int(line.rstrip("\n").lstrip("No "))

            #---------------------------------------------------------#
                    
            # If it is the first line of metadata about the alignment
            elif line.startswith(">"):
                    
                # Find the name of the hit 
                hit_name = line.split(" ")[0].lstrip(">")
                hit_name = "{0}.{1}".format(\
                           hit_name.split("_")[0].lower(),
                           hit_name.split("_")[1])
                    
                # Find the name of the organism with a regular
                # expression and get rid of the parentheses
                org = re.findall(r"{.*?}", \
                        line)[0].rstrip("}").lstrip("{").lower()

                # The first letter of the name of the organism
                # (first letter of the genus) must be uppercase
                org = org[0].upper() + org[1:]

            #---------------------------------------------------------#

            # If it is the second line of metadata about the alignment
            elif line.startswith("Probab"):
                    
                # Get all the data on it
                prob, evalue, score, aligned_cols, \
                identities, similarity, sum_probs, \
                template_neff = \
                    [float(item.split("=")[1].rstrip("%"))
                     for item in line.split(" ") if item != ""]
                    
                # Set the data to add to the corresponding hit in 
                # the table
                data_to_add = \
                    [identities, similarity, sum_probs,
                     template_neff, org]
                    
                # Generate a hit ID compatible with the one used
                # in the table
                hit_id = f"{hit_num}{self.HIT_NUM_ID_SEP}{hit_name}"
                    
                # Add the info from the alignment to the dictionary
                # storing them for all the hits
                info_from_als[hit_id] = \
                    dict(zip(self.COLNAMES_ALS, data_to_add))

                # Create the dictionary for the new alignment in the
                # 'als' dictionary
                als[hit_id] = \
                    {"Q": {}, "T": {}, 
                     "quality" : {},
                     "confidence" : {}}
                
                # Reset the sub-lines of the alignment
                sublines = []

            #---------------------------------------------------------#

            else:

                # Split the line by spaces and get rid of empty strings
                line = \
                    [item for item in line.rstrip("\n").split(" ") \
                     if item != ""]
                        
                # If the line starts with 'Q' or 'T', we are parsing
                # data about the query and the template
                if line[0] == "Q" or line[0] == "T":
                    
                    # - The section (e.g. 'Q', 'T') is the first item
                    # - The subsection (e.g. 'ss_pred', 'Consensus')
                    #   is the second item
                    # - The sequence is the third item
                    sec, subsec, seq = line[0:3]

                    #-------------------------------------------------#

                    # If the subsection is storing the aligned portion
                    # of the query
                    if subsec == "sp":
                        
                        # The subsection will be called 'query'
                        subsec = "query"

                    #-------------------------------------------------#
                    
                    # If the subsection is storing the aligned portion
                    # of the template
                    elif subsec[:4].lower() == hit_name[:4]:
                        
                        # The subsection will be called 'template'
                        subsec = "template"

                    #-------------------------------------------------#

                    # The third item is the residue counter in lines
                    # that have them
                    if seq.isdigit():

                        # If the section regards the query sequence
                        if sec == "Q":
                            
                            # If no record of the starting point for
                            # the query sequence was found so far (i.e.
                            # we are parsing the first portion of the
                            # alignment, where the start of the current
                            # portion is the start of the whole aligned
                            # portion)
                            if start_q is None:

                                # Save the starting point
                                start_q = int(line[2])

                            # Save the ending point of the current
                            # portion of the query sequence (it will
                            # be the ending point of the whole aligned
                            # portion only if the alignment is not
                            # split on multiple lines)  
                            end_q = int(line[4])

                            # If no info about the percent of the full
                            # UniProt query sequence covered by the
                            # current alignment was saved yet
                            if num_res is None:

                                # Get the total number of residues of
                                # the UniProt query sequence
                                num_res = \
                                    int(line[5].lstrip(\
                                        "(").rstrip(")"))

                        #---------------------------------------------#

                        # If the section regards the template sequence
                        elif sec == "T":
                            
                            # If no record of the starting point for
                            # the template sequence was found so far
                            # (i.e. we are parsing the first portion of
                            # the alignment, where the start of the
                            # current portion is the start of the whole
                            # aligned portion)
                            if start_t is None:

                                # Save the starting point
                                start_t = int(line[2])

                            # Save the ending point of the current
                            # portion of the template sequence (it
                            # will be the ending point of the whole
                            # aligned portion only if the alignment
                            # is not split on multiple lines)
                            end_t = int(line[4])

                        #---------------------------------------------#

                        # The sequence is the fourth item
                        seq = line[3]

                    #-------------------------------------------------#
                            
                    # If it is the first time we are parsing that
                    # sub-line (sub-lines are split across multiple
                    # chunks)
                    if not subsec in als[hit_id][sec].keys():

                        # Add a new string for the data
                        als[hit_id][sec][subsec] = seq

                    # Otherwise
                    else:

                        # Update the old one
                        als[hit_id][sec][subsec] += seq

                #-----------------------------------------------------#

                # If the line regards the confidence of the alignment
                elif line[0] == "Confidence":

                    # Set the name of the section
                    sec = "confidence"

                    # The sequence is the second item
                    seq = line[1]

                    # If it is the first time we are parsing the
                    # quality sub-line
                    if len(als[hit_id][sec]) == 0:

                        # Add a new string for it
                        als[hit_id][sec] = seq

                    # Otherwise
                    else:

                        # Update the old one
                        als[hit_id][sec] += seq

                #-----------------------------------------------------#

                # If the line does not start with either 'Q' or 'T' or
                # Confidence', we are parsing a line regarding the
                # alignment's quality
                else:

                    # Set the name of the section
                    sec = "quality"

                    # The sequence is the first item
                    seq = line[0]

                    # if it is the first time we are parsing the
                    # quality sub-line
                    if len(als[hit_id][sec]) == 0:

                        # Add a new string for it
                        als[hit_id][sec] = seq

                    # Otherwise
                    else:

                        # Update the old one
                        als[hit_id][sec] += seq
         
        #-------------------------------------------------------------#          

        # Convert the additional info found in the alignments to a
        # data frame
        df_als = pd.DataFrame.from_dict(info_from_als,
                                        orient = "index")

        #-------------------------------------------------------------#   

        # Return the alignments and the info about them
        return als, df_als


    def parse(self,
              hhsearch_output):
        """Parse the output file produced by HHsearch.

        Parameters
        ----------
        hhsearch_output : ``str``
            The output file produced by HHsearch.

        Returns
        -------
        df_final : ``pandas.DataFrame``
            A data frame containing information about the hits found
            by HHsearch.

        als : ``dict``
            A dictionary containing the alignments produced by
            HHsearch. 
        """

        # Open the HHsearch output file
        with open(hhsearch_output, "r") as f:

            # Get the indexes of the starting and ending lines of
            # the file's hits' and alignments' sections (it is an
            # open interval on the right side)
            range_hits, range_als = self._get_section_indexes(f)

            #---------------------------------------------------------#

            # Get the alignments and the info associated with them
            als, df_als = self._parse_als(f = f,
                                          range_indexes = range_als)

            #---------------------------------------------------------#

            # Get the hits' table
            df_hits = self._parse_hits(f = f,
                                       range_indexes = range_hits)

            #---------------------------------------------------------#

            # Concatenate the data frames containing information from
            # the hits and alignments sections
            df_final = pd.concat([df_hits, df_als],
                                 axis = 1,
                                 sort = False)

            #---------------------------------------------------------#

            # Return the final data frame and the alignments
            return df_final, als


#######################################################################


def get_aligned_region(hit_al,
                       seq_range):
    """Given a range of positions in the UniProt sequence, get the
    sequence's aligned portion in an alignment parsed from the output
    file produced by HHsearch.

    Parameters
    ----------
    hit_al : ``dict``
        One of the alignments parsed from the output file.

    seq_range : ``list`` or ``tuple``
        The starting and ending point of a sequence's portion to
        find in the alignment.

    Returns
    -------
    start : ``int``
        The starting point of the aligned sequence.

    seq : ``str``
        The aligned sequence.

    end : ``int``
        The ending point of the aligned sequence.

    is_match : ``bool``
        Whether the aligned sequence perfectly matches the
        original one.

    It returns ``None`` if the sequence's portion is not covered at
    all by the alignment.
    """

    # Get the starting end ending point of the sequence of interest
    start_s, end_s = seq_range

    #-----------------------------------------------------------------#
    
    # Get the sequence and the starting and ending points of the
    # aligned portion of the query 
    seq_q = hit_al["Q"]["query"]
    start_q = hit_al["Q"]["start"]
    end_q = hit_al["Q"]["end"]

    #-----------------------------------------------------------------#

    # Get the sequence and the starting and ending points of the
    # aligned portion of the template
    seq_t = hit_al["T"]["template"]
    start_t = hit_al["T"]["start"]
    end_t = hit_al["T"]["end"]

    #-----------------------------------------------------------------#
    
    # Get the full range of positions covered by the aligned portion
    # of the template 
    range_t = list(range(start_t, end_t))

    #-----------------------------------------------------------------#

    # If the sequence of interest is not covered by the aligned portion
    # of the query sequence
    if not (start_s >= start_q and end_s <= end_q):

        # Warn the user
        warnstr = \
            f"The sequence requested starts at position {start_s} " \
            f"and ends at position {end_s}, but the aligned portion " \
            "of the query in this alignment only covers positions " \
            f"{start_q} to {end_q}."
        logger.warning(warnstr)

        # Return
        return

    #-----------------------------------------------------------------#

    # Get the relative starting point of the sequence of interest
    # within the aligned portion of the query
    rel_start_s = start_s - start_q

    #-----------------------------------------------------------------#
    
    # Get the relative ending point of the sequence of interest
    # within the aligned portion of the query
    rel_end_s = rel_start_s + (end_s - start_s) + 1

    #-----------------------------------------------------------------#

    # Get the sequence in the template aligned to the sequence of
    # interest
    aligned_s = seq_t[rel_start_s:rel_end_s]

    #-----------------------------------------------------------------#

    # Get the starting and ending points of the aligned sequence in
    # the template
    start_s_in_t, end_s_in_t = \
        range_t[rel_start_s], range_t[rel_end_s]

    #-----------------------------------------------------------------#

    # If the aligned portion of the template matches the sequence
    # of interest
    if aligned_s == seq_q[rel_start_s:rel_end_s]:
        
        # Return the aligned portion of the template, its starting
        # and ending points, and 'True' to indicate that it matches
        # the sequence of interest
        return start_s_in_t, aligned_s, end_s_in_t, True

    #-----------------------------------------------------------------#

    # If the aligned portion of the template does not match the
    # sequence of interest
    else:
        
        # Return the aligned portion of the template, its starting
        # and ending points, and 'False' to indicate that it does not
        # match the sequence of interest
        return start_s_in_t, aligned_s, end_s_in_t, False
    