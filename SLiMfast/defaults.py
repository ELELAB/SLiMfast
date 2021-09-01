#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    defaults.py
#
#    Default values.
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
from pkg_resources import resource_filename, Requirement


# Metadata (taken from MMCIF files) retrieved for PDB IDs
# of interest
STRUCTURE_METADATA = [
    # Experimental method used for structure determination
    ("_exptl.method",
     str, 
     "experimental_method"),
    # Resolution
    ("_refine.ls_d_res_high",
     float,
     "resolution"),
    # Deposition date
    ("_pdbx_database_status.recvd_initial_deposition_date",
     str,
     "deposition_date"),
    ]

# Separator for starting and ending points of different
# UniProt segments
SEGMENTS_SEP = ";"

# Separator for numeric ranges
RANGE_SEP = "-"

# Data frame field containing PDB IDs
PDB_ID_FIELD = "pdb_id"

# Data frame field containing PDB chains
PDB_CHAIN_FIELD = "pdb_chain"

# Data frame field containing UniProt IDs
UNP_ID_FIELD = "uniprot_id"

# Data frame field containing ranges representing
# UniProt segments 
UNP_SEGMENTS_FIELD = "uniprot_segments"

# Data frame field containing ranges representing
# segments of PDB structures
PDB_SEGMENTS_FIELD = "pdb_segments"