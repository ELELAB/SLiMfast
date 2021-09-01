# SLiMfast

**Last updated**: September 1st, 2021

## Introduction

`SLiMfast` is a Python package implementing a simple pipeline to annotate a set of putative Short Linear Motifs (SLiMs) in a protein sequence according to whether they fall in structured or unstructured regions, using structural information available in the Protein Data Bank (PDB).

If no structures of the protein of interest are found (or they do not cover the SLiMs under consideration), `SLiMfast` provides a tool to search for homologous proteins with available structures in the PDB to be able to annotate the SLiMs. To do so, `SLiMfast` relies on the `hhblits` and `hhsearch` tools of the `HHsuite`, which must be installed in order for `SLiMfast` to work.

In case no structures of either the protein of interest or of any homologs are available in the PDB, `SLiMfast` will interrogate the `DisProt` database to find annotations of disorder for the SLiMs of interest.

Since SLiMs often falls within unstructured regions, `SLiMfast` provides insights about which of the SLiMs are likely to be functional in the context of the whole protein.

## Installation

Detailed instructions to install `SLiMfast` system-wide or in virtual environments (with `virtualenv` or `conda`) are provided in the `INSTALL.md` file.