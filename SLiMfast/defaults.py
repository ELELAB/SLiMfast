#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    defaults.py
#
#    Default values.
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
from pkg_resources import resource_filename, Requirement



# Data frame field containing a structure's AlphaFold ID
ALPHAFOLD_ID_FIELD = "alphafold_id"

# Directory where the substitution matrices are stored
CONFIG_RUN_DIR = \
    resource_filename(Requirement("SLiMfast"), "SLiMfast/config_run")

# Data frame field containing the unique hit identifier
# for the homologs
HIT_ID_FIELD = "hit_id"

# Separator for numeric ranges
RANGE_SEP = "-"

# Data frame field containing the deposition date
DEPOSITION_DATE_FIELD = "deposition_date"

# Data frame field containing the method used to obtain
# the structure
METHOD_FIELD = "method"

# Data frame field containing the resolution of the
# structure (if applicable)
RESOLUTION_FIELD = "resolution"

# Separator for starting and ending points of different
# UniProt segments
SEGMENTS_SEP = ";"

# Data frame field containing the structures' chains
STRUCT_CHAIN_FIELD = "struct_chain"

# Data frame field containing the fragment number of
# a structure from the AlphaFold Protein Structure
# Database
STRUCT_FRAGMENT_FIELD = "struct_fragment"

# Data frame field containing the structures' IDs
STRUCT_ID_FIELD = "struct_id"

# Data frame field containing ranges representing
# segments of the structures
STRUCT_SEGMENTS_FIELD = "struct_segments"

# Data frame field containing the version of
# a structure from the AlphaFold Protein Structure
# Database
STRUCT_VERSION_FIELD = "struct_version" 

# Data frame field containing UniProt IDs
UNP_ID_FIELD = "uniprot_id"

# data frame field containing the UniProt isoform
UNP_ISOFORM_FIELD = "uniprot_isoform"

# Data frame field containing ranges representing
# UniProt segments 
UNP_SEGMENTS_FIELD = "uniprot_segments"

