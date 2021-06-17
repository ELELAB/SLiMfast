#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    run_tools.py
#
#    Utilities to run external tools needed for the SLiMfast pipeline.
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
import os
import subprocess
# SLiMfast
from SLiMfast import (
    uniprot_tools,
    util,
    )


def run_process(executable,
                options, 
                wd, 
                dirname, 
                outfile = None):
    """Run a generic process."""

    # Make sure that the directory where to run the process exists.
    # If not, create it.
    os.makedirs(os.path.join(wd, dirname), exist_ok = True)

    # The output of the process will be written either to
    # the standard output or to an output file
    out = \
        open(outfile, "w") if outfile is not None \
        else subprocess.STDOUT

    # Launch the process (redirect errors to the standard output)
    popen = subprocess.Popen([executable] + args,
                             stdout = out,
                             stderr = subprocess.STDOUT,
                             cwd = wd)
    
    # Wait for the process to complete
    popen.wait()
    
    # Return the Popen attributes of interest (cannot
    # return Popen itself since it is not serializable and
    # we are launching the process with Dask)
    return {"args" : popen.args, \
            "stdin" : popen.stdin, \
            "stdout" : popen.stdout, \
            "stderr" : popen.stderr, \
            "pid" : popen.pid, \
            "returncode" : popen.returncode}


def run_hhpred(uniprot_id,
               hhblits_params,
               hhsearch_params,
               wd):
    """Run the HHpred pipeline to find remote homologs of a
    a protein from an input multiple sequence alignment (MSA)
    built from the protein sequence.
    """


    #----------------------- Get UniProt FASTA -----------------------#


    # Set the path for the FASTA file of the sequence to be written
    uniprot_fasta = os.path.join(wd, uniprot_id + ".fasta")
    
    # Get the sequence corresponding to the UniProt ID and save
    # it to a FASTA file
    uniprot_tools.write_uniprot_fasta(uniprot_id = uniprot_id,
                                      fasta_file = uniprot_fasta)

    # Get header and sequence from the FASTA file
    header, sequence = util.read_fasta(fasta_file = uniprot_fasta)


    #-------------------------- Run HHblits --------------------------#


    # Make a copy of the parameters set to run HHblits
    # (so that the original dictionary is not modified)
    hhblits_params_copy = dict(hhblits_params)

    # Path to the directory where HHblits will be run
    hhblits_dirname = os.path.join(wd, hhblits_params["dirname"])

    # Set the path to the directory where HHblits will be run
    hhblits_params_copy["dirname"] = hhblits_dirname

    # Set the path to the output A3M file
    uniprot_a3m = \
        os.path.join(hhblits_dirname, uniprot_id + ".a3m")
            
    # Add the FASTA file as input and the A3M file as output
    # to the options of HHblits
    hhblits_params_copy["options"] = \
        {**hhblits_params["options"],
         "-i" : uniprot_fasta, "-o3am" : uniprot_a3m}
        
    # Run HHblits for the generation of the query MSA
    hhblits_process = run_process(**hhblits_params_copy)


    #-------------------------- Run HHsearch -------------------------#


    # Make a copy of the parameters set to run HHsearch
    hhsearch_params_copy = dict(hhsearch_params)

    # Path to the directory where HHsearch will be run
    hhsearch_dirname = os.path.join(wd, hhsearch_params["dirname"])

    # Set the path to the directory where HHsearch will be run
    hhsearch_params_copy["dirname"] = hhsearch_dirname

    # Path to the HHsearch output
    uniprot_hhsearch = \
        os.path.join(hhsearch_dirname, uniprot_id + ".hhsearch")

    # Add the A3M file as input
    hhsearch_params_copy["options"] = \
        {**hhsearch_params_copy["options"],
         **{"-i" : uniprot_a3m, "-o" : uniprot_hhsearch}}
        
    # Search profile HHM database(s) with HHsearch
    hhsearch_process = run_process(**hhsearch_params_copy)

    return (hhblits_process, hhsearch_process), hhsearch_output


def run_dssp(pdb_file, dssp_params, wd):
    """Run DSSP to extract secondary structure information
    from a PDB file.
    """

    # Make a copy of the parameters set to run DSSP
    dssp_params_copy = dict(dssp_params)

    # Set the path to the DSSP output file
    dssp_output = \
        os.path.join(wd,
                     os.path.basename(pdb_file)[:-4] + ".dssp")
                
    # Add the PDB file as input to the options for DSSP
    dssp_params_copy["options"] = \
        {**dssp_params_copy["options"],
         **{"-i" : pdb_file, "-o" : dssp_output}}
                
    # Run DSSP
    dssp_process = run_process(**dssp_params_copy)

    return dssp_process, dssp_output