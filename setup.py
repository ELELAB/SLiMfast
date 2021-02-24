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
package_data = {"SLiMfast" : ["config_plot/*",
                              "config_run/*",
                              "config_settings/*"]}
entry_points = \
    {"console_scripts" : \
        ["slimfast_run = SLiMfast.slimfast_run:main", \
         "slimfast_plot = SLiMfast.slimfast_plot:main"], \
    }
install_requires = ["biopython", \
                    "dask", \
                    "distributed", \
                    "matplotlib", \
                    "MDAnalysis", \
                    "numpy", \
                    "pandas", \
                    "pyyaml", \
                    "seaborn"]

setup(name = name, \
      url = url, \
      author = author, \
      version = version, \
      description = description, \
      include_package_data = True, \
      package_data = package_data, \
      package_dir = package_dir, \
      packages = packages, \
      entry_points = entry_points, \
      install_requires = install_requires)