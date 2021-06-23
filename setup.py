#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

from setuptools import setup

name = "SLiMfast"
url = "https://github.com/ELELAB/SLiMfast"
author = "Valentina Sora, Matteo Tiberti, Elena Papaleo"
version = "0.0.1"
description = "Structural filtering of putative SLiMs in a protein sequence."
package_dir = {"SLiMfast" : "SLiMfast"}
packages = ["SLiMfast"]
entry_points = \
    {"console_scripts" : \
        ["slimfast_get_structures = SLiMfast.slimfast_get_structures:main",
         "slimfast_get_homologs = SLiMfast.slimfast_get_homologs:main",
         "slimfast_find_slims = SLiMfast.slimfast_find_slims:main"], \
    }
install_requires = ["biopython",
                    "pandas",
                    "pytest",
                    "pyyaml",
                    "requests"]

setup(name = name,
      url = url,
      author = author,
      version = version,
      description = description,
      include_package_data = True,
      package_dir = package_dir,
      packages = packages,
      entry_points = entry_points,
      install_requires = install_requires)