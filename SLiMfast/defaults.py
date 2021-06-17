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


STRUCTURE_METADATA = [
    # experimental method used for structure determination
    ("_exptl.method",
     str, 
     "experimental_method"),
    # resolution
    ("_refine.ls_d_res_high",
     float,
     "resolution"),
    # deposition date
    ("_pdbx_database_status.recvd_initial_deposition_date",
     str,
     "deposition_date"),
    ]

# Separator for starting and ending points of different
# UniProt segments
SEGMENTS_SEP = ";"

RANGE_SEP = "-"

PDB_ID_FIELD = "pdb_id"
PDB_CHAIN_FIELD = "pdb_chain"
UNP_ID_FIELD = "uniprot_id"
UNP_SEGMENTS_FIELD = "uniprot_segments"
PDB_SEGMENTS_FIELD = "pdb_segments"