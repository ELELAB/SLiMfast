#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    setup.py
#
#    SLiMfast setup.
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
from setuptools import setup


# Package name
name = "SLiMfast"

# GitHub URL
url = "https://github.com/ELELAB/SLiMfast"

# Author(s)
author = "Valentina Sora, Matteo Tiberti, Elena Papaleo"

# Package version
version = "0.0.1"

# Short description of the package
description = \
    "Structural filtering of putative SLiMs in a protein sequence."

# Mapping between the package name and the corresponding directoru
package_dir = {"SLiMfast" : "SLiMfast"}

# Packages included
packages = ["SLiMfast"]

# Package data to include
package_data = {"SLiMfast" : ["config_run/*"]}

# Entry points
entry_points = \
    {"console_scripts" : \
        ["slimfast_get_structures = SLiMfast.slimfast_get_structures:main",
         "slimfast_get_homologs = SLiMfast.slimfast_get_homologs:main",
         "slimfast_find_slims = SLiMfast.slimfast_find_slims:main"], \
    }

# Dependencies required
install_requires = ["biopython",
                    "pandas",
                    "pytest",
                    "pyyaml",
                    "requests"]

# Run the setup
setup(name = name,
      url = url,
      author = author,
      version = version,
      description = description,
      package_dir = package_dir,
      packages = packages,
      include_package_data = True,
      package_data = package_data,
      entry_points = entry_points,
      install_requires = install_requires)