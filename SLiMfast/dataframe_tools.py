#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    dataframe_tools.py
#
#    Utility functions to manipulate Pandas data frames.
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



def filter_df(df,
              filters):
    """Filter the data frame.
    """

    # For each filter passed
    for arg, val in filters.items():

        # If it defines a minimum value
        if arg.startswith("min_"):
            
            # The name of the column it refers to will be
            # after min_
            col = arg[4:]
            
            # Set the filtering condition
            cond = (df[col] >= val)
                
        # If it defines a maximum value
        elif arg.startswith("max_"):
            
            # The name of the column it refers to will be
            # after max_
            col = arg[4:]
            
            # Set the filtering condition
            cond = (df[col] <= val)

        # If it defines a series of values that can be
        # matches
        elif arg.endswith("_in"):

            # The name of the column it refers to will be
            # before _in
            col = arg[:-3]

            # Set the filtering condition
            cond = df[col].isin(val)

        # If it defines an exact value (the argument is the
        # name of a column)
        elif arg in df.columns:
            
            # Set the filtering condition
            cond = (df[arg] == val)

        # Otherwise, it is an invalid keyword
        else:

            # Inform the user and raise an error
            errstr = \
                f"Invalid keyword '{arg}' passed for filtering " \
                f"the data frame."
            raise ValueError(errstr)

        # Filter the data frame
        df = df[cond]

    # Return the filtered data frame
    return df