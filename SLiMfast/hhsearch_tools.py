#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    hhseach_tools.py
#
#    Contains a class implementing a parser to parse data from 
#    HHsearch output files and utility functions to manipulate
#    such data.
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
import re
# Third-party packages
import pandas as pd



############################### CLASSES ###############################



class HHsearchParser:

    """Class implementing a parser for HHsearch output files.
    """

    # Format for data in scientific notation
    SCI_FMT = "{:,.1E}"

    # Column names for the hits' table
    COLNAMES_HITS = ["num", "hit", "description", "prob",
                     "e_value", "p_value", "score", "ss", "cols",
                      "uniprot_segments", "pdb_segments"]

    # Column names for the data retrieved from the aligments
    COLNAMES_ALS = ["identities", "similarity", "sum_probs",
                    "template_neff", "organism"]

    # Column names for the additional columns that will
    # contain PDB ID and chain
    COLNAMES_PDB = ["pdb_id", "pdb_chain"]

    # Column name for the additional column that will
    # contain information about the sequence coverage
    COLNAME_COV = "uniprot_coverage_percent"

    # Separator to be used between the hit number and the
    # hit ID (number + ID will be the hit's unique
    # identifier)
    HIT_NUM_ID_SEP = "_"


    def __init__(self):
        """Initialize the parser.
        """
        pass


    def _get_section_indexes(self, f):
        """Get the indexes of the hits' and alignments' sections.
        """

        # Boolean switch to indicate when to parse the table
        # containing the hits; off at the beginning since
        # the table comes after the metadata
        parse_hits = False
            
        # Boolean switch to indicate when to parse the alignments;
        # off at the beginning since they come after the hits' table
        parse_als = False

        # Initialize the indexes of the lines where the hits'
        # table starts and ends
        start_hits_ix = 0
        end_hits_ix = 0
        
        # For each line
        for ix, line in enumerate(f):

            # If you encounter an empty line
            if re.match(r"^\s*$", line):

                if not parse_hits and not parse_als:
                    
                    # Turn on the parsing for the hits' table
                    # since it is the next section
                    parse_hits = True
                    
                    # The starting index for the hits' table
                    # section is the one of the next line
                    start_hits_ix = ix + 1
                    
                    # Ignore the current (empty) line
                    continue

            # If you are parsing the hits' table section
            if parse_hits:

                # If you encounter a blank line
                if re.match(r"^\s*$", line):
                    
                    # Turn off the parsing of the hits' table
                    parse_hits = False
                    
                    # Turn on the parsing for the alignments
                    parse_als = True
                    
                    # The ending index for the hits' table
                    # is the one of the current line
                    # (the interval is open, so the table
                    # ends on the previous line)
                    end_hits_ix = ix
                    
                    # Ignore the current (empty) line
                    continue

        # Returns the indexes of the hits' and the aligments'
        # sections
        return (start_hits_ix, end_hits_ix), (end_hits_ix+1, ix)


    def _parse_hits(self, f, range_indexes):
        """Parse the hits' section.
        """

        # Create an empty dictionary to store data
        data = []

        # Boolean switch for when to start parsing the hits'
        # section 
        parse = False

        # Unpack the range in starting and ending index
        start_ix, end_ix = range_indexes

        # Go to the beginning of the file
        f.seek(0)
        
        # For each line
        for i, line in enumerate(f):

            # If you reached the beginning of the hits' section,
            # start parsing
            if i == start_ix + 1:
                parse = True

            # If you reached the end of the hits' section,
            # exit from the loop
            elif i == end_ix:
                break

            # If parsing is not active, continue
            if not parse:
                continue

            # Split the line with whitespaces and take all
            # non empty elements
            line = [i for i in line.split(" ") if i != ""]

            # Get the hit ID
            hit = \
                line[1].split("_")[0].lower() + "." \
                + line[1].split("_")[1]

            # If the PDB segments are stored in the last
            # elements of the line
            if "-" in line[-1]:
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
                     "uniprot_segments" : line[-2],
                     "pdb_segments" : line[-1].split("(")[0]})

            # If the PDB segments are stored in the second-to-last
            # element of the line
            else:
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
                     "uniprot_segments" : line[-3],
                     "pdb_segments" : line[-2]})

        # Convert the dictionary containing the hits into a data frame
        hits_df = pd.DataFrame(data)

        # Specify the column names
        hits_df.columns = self.COLNAMES_HITS

        # Add two columns storing just the PDB ID and the chain
        hits_df[self.COLNAMES_PDB] = \
            hits_df["hit"].str.split(".", 1, expand = True)

        # Row names are composed by both the hit number and the
        # PDB ID (with chain) to both ensure their uniqueness and
        # be informative at the same time
        hits_df.index = [str(n) + self.HIT_NUM_ID_SEP + h for n, h \
                         in zip(hits_df["num"], hits_df["hit"])]

        # Return the data frame
        return hits_df


    def _parse_als(self, f, range_indexes):
        """Parse the alignment sections and returns a dictionary of
        dictionaries containing the strings corresponding to
        each alignment. Having the following alignment in a file
        the dictionary will be constructed as shown below:
        
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

        # Set the boolean switch that indicates to
        # parse the alignment to False
        parse = False

        # Set the starting and ending points of query
        # and template to None
        start_q = None
        start_t = None
        end_q = None
        end_t = None

        # Set the hit ID to None
        hit_id = None

        # Set the number of residues in the query to None
        num_res = None

        # Initialize an empty dictionary to store the info
        # on the hits gathered from the alignments
        info_from_als = {}

        # Create an empty dictionary to atore the alignments
        als = {}
        
        # Go to the beginning of the file
        f.seek(0)

        # For each line
        for ix, line in enumerate(f):

            # If you reached the beginning of the alignments'
            # section, start parsing
            if ix == range_indexes[0]:
                parse = True

            # If parsing is not active, continue
            if not parse:
                continue

            # Ignore empty lines
            if re.match(r"^\s*$", line):
                continue

            # If you just entered a new aligment...
            elif line.startswith("No"):

                # If there are data from a previous hit to
                # be stored
                if hit_id is not None:

                    # If no info for the coverage of the UniProt
                    # sequence by the alignment has been stored
                    # for the hit
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
                    
                    # Save the start and end numbers of both
                    # the query and the template/target
                    als[hit_id]["Q"]["start"] = start_q
                    als[hit_id]["Q"]["end"] = end_q
                    als[hit_id]["T"]["start"] = start_t
                    als[hit_id]["T"]["end"] = end_t

                    # Re-initialize all of them to None for
                    # the next hit
                    start_q = None
                    end_q = None
                    start_t = None
                    end_t = None
                
                # Get the corresponding HHpred hit number
                hit_num = int(line.rstrip("\n").lstrip("No "))

                # Go to the next line
                continue
                    
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

                continue

            # If it is the second line of metadata about the alignment
            elif line.startswith("Probab"):
                    
                # Get all the data on it
                prob, evalue, score, aligned_cols, \
                identities, similarity, sum_probs, \
                template_neff = \
                    [float(item.split("=")[1].rstrip("%"))
                     for item in line.split(" ") if item != ""]
                    
                # Data to add to the corresponding hit in the table
                data_to_add = [identities, similarity, sum_probs,
                               template_neff, org]
                    
                # Generate a hit ID compatible with the one used
                # in the table
                hit_id = f"{hit_num}{self.HIT_NUM_ID_SEP}{hit_name}"
                    
                # Add the info from the alignment to the dictionary
                # storing them for all the hits
                info_from_als[hit_id] = \
                    dict(zip(self.COLNAMES_ALS, data_to_add))

                # Create the dictionary for the new 
                # alignment in the 'als' dictionary
                als[hit_id] = {"Q": {}, "T": {}, 
                               "quality" : {}, "confidence" : {}}
                
                # Reset the sub-lines of the alignment
                sublines = []

                # Continue to the next line
                continue

            # Split the line by spaces and get rid of empty strings
            line = [item for item in line.rstrip("\n").split(" ") \
                    if item != ""]
                    
            # If the line starts with 'Q' or 'T', we are
            # parsing data about the query and the template
            if line[0] == "Q" or line[0] == "T":
                
                # The section (e.g. 'Q', 'T') is the first item

                # The subsection (e.g. 'ss_pred', 'Consensus')
                # is the second item

                # The sequence is the third item

                sec, subsec, seq = line[0:3]

                # If the subsection is storing the aligned
                # portion of the query
                if subsec == "sp":
                    # The subsection will be called 'query'
                    subsec = "query"
                # If the subsection is storing the aligned
                # portion of the template
                elif subsec[:4].lower() == hit_name[:4]:
                    # The subsection will be called 'template'
                    subsec = "template"

                # The third item is the residue counter in lines
                # that have them
                if seq.isdigit():

                    # If the section regards the query sequence
                    if sec == "Q":
                        # If no record of starting point for the
                        # query sequence was found so far (i.e.
                        # we are parsing the first portion of
                        # the alignment, where the start of the
                        # current portion is the start of the
                        # whole aligned portion), save the
                        # starting point
                        if start_q is None:
                            start_q = int(line[2])

                        # Save the ending point of the current
                        # portion of the query sequence (will
                        # be the ending point of the whole aligned
                        # portion only if the alignment is not
                        # split on multiple lines)  
                        end_q = int(line[4])

                        # If no info about the percent of the full UniProt
                        # query sequence covered by the current alignment
                        # has been saved yet

                        if num_res is None:
                            # Get the total number of residues of the
                            # UniProt query sequence
                            num_res = int(line[5].lstrip("(").rstrip(")"))

                    # If the section regards the template sequence
                    elif sec == "T":
                        # If no record of starting point for the
                        # template sequence was found so far (i.e.
                        # we are parsing the first portion of
                        # the alignment, where the start of the
                        # current portion is the start of the
                        # whole aligned portion), save the
                        # starting point
                        if start_t is None:
                            start_t = int(line[2])

                        # Save the ending point of the current
                        # portion of the template sequence (will
                        # be the ending point of the whole aligned
                        # portion only if the alignment is not
                        # split on multiple lines)
                        end_t = int(line[4])

                    # The sequence is the fourth item
                    seq = line[3]

                        
                # If it is the first time we are parsing that
                # sub-line (sub-lines are splitted across
                # multiple chuncks), add a new string for
                # the data
                if not subsec in als[hit_id][sec].keys():
                    als[hit_id][sec][subsec] = seq

                # ...otherwise, update the old one
                else:
                    als[hit_id][sec][subsec] += seq


            # If the line regards the confidence of the alignment
            elif line[0] == "Confidence":

                # Name of the section
                sec = "confidence"

                # The sequence is the second item
                seq = line[1]

                # if it is the first time we are parsing the
                # quality sub-line, add a new string for it
                if len(als[hit_id][sec]) == 0:
                    als[hit_id][sec] = seq

                # ... otherwise, update the old one
                else:
                    als[hit_id][sec] += seq


            # If the line does not start with either 'Q' or 
            # 'T' or 'Confidence', we are parsing a line
            # sorting the alignment quality
            else:

                # Name of the section
                sec = "quality"

                # The sequence is the first item
                seq = line[0]

                # if it is the first time we are parsing the
                # quality sub-line, add a new string for it
                if len(als[hit_id][sec]) == 0:
                    als[hit_id][sec] = seq

                # ... otherwise, update the old one
                else:
                    als[hit_id][sec] += seq

        # Convert the additional info found in the alignments
        # to a dataframe
        als_df = pd.DataFrame.from_dict(info_from_als,
                                        orient = "index")

        # Return the alignments and the info about them
        return als, als_df


    def parse(self, hhsearch_output):
        """Parse the HHpred output file.
        """

        with open(hhsearch_output, "r") as f:

            # Get the indexes of the starting end ending
            # line of its section (open interval on the
            # right side)
            range_hits, range_als = self._get_section_indexes(f)

            # Get the alignments and the info associated with them
            als, als_df = self._parse_als(f = f,
                                          range_indexes = range_als)

            # Get the hits' table
            hits_df = self._parse_hits(f = f,
                                       range_indexes = range_hits)

            # Concatenate the dataframes containing information from
            # the hits and alignments sections
            final_df = pd.concat([hits_df, als_df],
                                 axis = 1,
                                 sort = False)

            # Return the final dataframe and the alignments themselves
            return final_df, als



########################## UTILITY FUNCTIONS ##########################



def filter_hits(hits, filters):
    """Filter the hits' dataframe.
    """

    # For each filter passed
    for arg, val in filters.items():
                
        # If it defines a minimum value
        if arg.startswith("min_"):
            # The name of the column it refers to will be
            # after min_
            col = arg[4:]
            # Set the filtering condition
            cond = (hits[col] >= val)
                
        # If it defines a maximum value
        elif args.startswith("max_"):
            # The name of the column it refers to will be
            # after max_
            col = arg[4:]
            # Set the filtering condition
            cond = (hits[col] <= val)

        # If it defines an exaxt value (the argument is the
        # name of a column)
        elif arg in cols:
            # Set the filtering condition
            cond = (hits[arg] == val)

        # Filter the dataframe
        hits = hits[cond]

    # Return the filtered dataframe
    return hits


def sort_hits(hits, sort_by, ascending):
    """Sort the hits' dataframe.
    """

    # Default to ascending sorting if no order is specified
    ascending = ascending if ascending is not None \
                else [True]*len(sort_by)
            
    # Sort the dataframe
    hits = hits.sort_values(by = sort_by,
                            ascending = ascending)

    # Return the sorted dataframe
    return hits


def get_aligned_region(hit_al, seq_range):
    """Given a range of positions in the UniProt sequence,
    get the aligned portion in the alignment.
    """

    # Get the starting end ending point of the sequence
    # of interest
    start_s, end_s = seq_range
    
    # Get the sequence and the starting and ending points
    # of the aligned portion of the query 
    seq_q = hit_al["Q"]["query"]
    start_q = hit_al["Q"]["start"]
    end_q = hit_al["Q"]["end"]

    # Get the sequence and the starting and ending points
    # of the aligned portion of the template
    seq_t = hit_al["T"]["template"]
    start_t = hit_al["T"]["start"]
    end_t = hit_al["T"]["end"]
    
    # Get the full range of positions covered by the
    # aligned portion of the template 
    range_t = list(range(start_t,end_t))

    # If the sequence of interest is not covered by the
    # aligned portion of the query sequence, return None
    if not (start_s >= start_q and end_s <= end_q):
        warnstr = \
            f"The sequence you requested starts at position " \
            f"{start_s} and end at position {end_s}, but the " \
            f"aligned portion of the query in this alignment " \
            f"only covers the positions {start_q}-{end_q}."
        log.warning(warnstr)
        return None

    # Get the relative starting point of the sequence of
    # interest within the aligned portion of the query
    rel_start_s = start_s - start_q
    
    # Get the relative ending point of the sequence of
    # interest within the aligned portion of the query
    rel_end_s = rel_start_s + (end_s - start_s) + 1

    # Get the sequence in the template aligned to
    # the sequence of interest
    aligned_s = seq_t[rel_start_s:rel_end_s]

    # Get the starting and ending points of the aligned
    # sequence in the template
    start_s_in_t, end_s_in_t = \
        range_t[rel_start_s], range_t[rel_end_s]

    # If the aligned portion of the template is identical
    # to the sequence of interest
    if aligned_s == seq_q[rel_start_s:rel_end_s]:
        # Return the aligned portion of the template, its
        # starting and ending points and 'True' to indicate
        # it perfectly matches the sequence of interest
        return start_s_in_t, aligned_s, end_s_in_t, True

    # If the aligned portion of the template does not match
    # the sequence of interest
    else:
        # Return the aligned portion of the template, its
        # starting and ending points and 'False' to indicate
        # it does not match the sequence of interest
        return start_s_in_t, aligned_s, end_s_in_t, False
    
