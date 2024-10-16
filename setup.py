#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    setup.py
#
#    SLiMfast setup.
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


# Standard library
from setuptools import setup


#######################################################################


# The name of the package
name = "SLiMfast"

# The URL where to find the package
url = f"https://github.com/ELELAB/{name}"

# The author(s) of the package
author = "Valentina Sora, Matteo Tiberti, Elena Papaleo"

# The version of the package
version = "0.0.2"

# A short description of the package
description = \
    "Structural filtering of putative SLiMs in protein sequences."

# A mapping between the package names and the corresponding directories
package_dir = {name : name}

# Which packages are included
packages = [name]

# Which package data to include
package_data = {name : ["config_run/*"]}

# The entry points
entry_points = \
    {"console_scripts" : \
        ["slimfast_get_structures = " \
            f"{name}.slimfast_get_structures:main",
         "slimfast_get_alphafold = " \
            f"{name}.slimfast_get_alphafold:main",
         "slimfast_get_homologs = " \
            f"{name}.slimfast_get_homologs:main",
         "slimfast_find_slims = " \
            f"{name}.slimfast_find_slims:main",
         "slimfast_run_pipeline = " \
            f"{name}.slimfast_run_pipeline:main"], \
    }

# The dependencies required
install_requires = ["biopython",
                    "dask",
                    "distributed",
                    "pandas",
                    "pytest",
                    "pyyaml",
                    "requests"]


#######################################################################


# Run the setup
setup(name = name,
      url = url,
      author = author,
      version = version,
      description = description,
      package_dir = package_dir,
      packages = packages,
      package_data = package_data,
      entry_points = entry_points,
      install_requires = install_requires)
