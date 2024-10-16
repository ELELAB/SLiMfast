# SLiMfast

Last updated: 17/04/2024

## Introduction

SLiMfast is a Python package implementing a simple pipeline to annotate a set of putative Short Linear Motifs (SLiMs) in a protein sequence according to whether they fall in structured or unstructured regions, using structural information available in the Protein Data Bank (PDB) and in the AlphaFold Protein Structure Database.

Suppose no structures of the protein of interest are found (or they do not cover the SLiMs under consideration). In that case, SLiMfast provides a tool to search for homologous proteins with available structures in the PDB to be able to annotate the SLiMs. To do so, SLiMfast relies on the HHblits and HHsearch tools of the HHsuite, which must be installed in order for SLiMfast to work.

If no structures of either the protein of interest or of any homologs are available, SLiMfast will interrogate the DisProt database to find disorder annotations for the SLiMs of interest.

Since SLiMs often fall within unstructured regions, SLiMfast provides insights about which SLiMs are likely to be functional in the context of the whole protein by using the DSSP program. SLiMfast can also calculate solvent accessibilities for specific residues using NACCESS.

The PDF user guide provided with the package provides a complete description of SLiMfast's functionalities and examples of how to use it.

## Installation

The `INSTALL.md` file provides detailed instructions for installing SLiMfast system-wide or in virtual environments (with `virtualenv` or `conda`).