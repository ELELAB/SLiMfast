#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    run_tools.py
#
#    Utilities to run external tools needed for the SLiMfast pipeline.
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
import logging as log
import os
import subprocess
# SLiMfast
from SLiMfast import (
    uniprot_tools,
    util,
    )



# Get the module logger
logger = log.getLogger(__name__)



def run_process(executable,
                options, 
                wd, 
                outfile = None):
    """Run a generic process."""

    # Make sure that the directory where to run the process exists.
    # If not, create it.
    os.makedirs(wd, exist_ok = True)

    # Convert the dictionary of options into a list of arguments
    args = \
        [str(item) for opt_val_pair in options.items() \
         for item in opt_val_pair]

    # If no output file has been specified
    if outfile == "" or outfile is None:
        
        # Launch the process (the output will be written to the
        # standard output, and errors will be redirected
        # to the standard output as well)
        popen = subprocess.Popen([executable] + args,
                                 stdout = subprocess.PIPE,
                                 cwd = wd)
    
    # If an output file has been specified
    else:

        # Open the file
        outfile = open(outfile, "w")

        # Launch the process (the output will be written to the
        # output file, and errors will be redirected to the
        # output file as well)
        popen = subprocess.Popen([executable] + args,
                                 stdout = outfile,
                                 stderr = outfile,
                                 cwd = wd)    
    
    # Wait for the process to complete
    popen.wait()
    
    # Return the Popen attributes of interest (cannot
    # return Popen itself since it is not serializable and
    # we may launch the process with Dask in the future)
    return {"args" : popen.args,
            "stdin" : popen.stdin,
            "stdout" : popen.stdout,
            "stderr" : popen.stderr,
            "pid" : popen.pid,
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
    uniprot_tools.write_uniprot_fasta(\
        uniprot_id = uniprot_id,
        fasta_file_handle = open(uniprot_fasta, "w"))

    # Get header and sequence from the FASTA file
    header, sequence = util.read_fasta(fasta_file = uniprot_fasta)


    #-------------------------- Run HHblits --------------------------#


    # Make a copy of the parameters set to run HHblits
    # (so that the original dictionary is not modified)
    hhblits_params_copy = dict(hhblits_params)

    # Set the path to the directory where HHblits will be run
    hhblits_dirname = os.path.join(wd, hhblits_params["wd"])

    # Set the path to the directory where HHblits will be run
    hhblits_params_copy["wd"] = hhblits_dirname

    # Set the path to the output A3M file
    uniprot_a3m = \
        os.path.join(hhblits_dirname, uniprot_id + ".a3m")

    # Add the FASTA file as input and the A3M file as output
    # to the options of HHblits
    hhblits_params_copy["options"] = \
        {**hhblits_params["options"],
         "-i" : uniprot_fasta, "-oa3m" : uniprot_a3m}

    # If the output file does not exist yet (no previous run
    # has been performed)
    if not os.path.exists(uniprot_a3m):

        # Inform the user that HHblits is about to be run
        infostr = \
            f"Running HHblits for the protein with UniProt ID " \
            f"{uniprot_id}..."
        logger.info(infostr)
        
        # Run HHblits for the generation of the query MSA
        hhblits_process = run_process(**hhblits_params_copy)
    
    # If the output file exists, do not run HHblits
    else:

        # Inform the user
        warnstr = \
            f"The file {uniprot_a3m} was found in {hhblits_dirname}." \
            f" SLiMfast will use it as input to HHsearch, without " \
            f"running HHblits again."
        logger.warning(warnstr)

        # No data will be associated with the HHblits process in this
        # case, since HHblits was not run
        hhblits_process = None


    #-------------------------- Run HHsearch -------------------------#


    # Make a copy of the parameters set to run HHsearch
    hhsearch_params_copy = dict(hhsearch_params)

    # Set the path to the directory where HHsearch will be run
    hhsearch_dirname = os.path.join(wd, hhsearch_params["wd"])

    # Set the path to the directory where HHsearch will be run
    hhsearch_params_copy["wd"] = hhsearch_dirname

    # Path to the HHsearch output
    uniprot_hhsearch = \
        os.path.join(hhsearch_dirname, uniprot_id + ".hhsearch")

    # Add the A3M file as input
    hhsearch_params_copy["options"] = \
        {**hhsearch_params_copy["options"],
         **{"-i" : uniprot_a3m, "-o" : uniprot_hhsearch}}

    # If the output file does not exist yet (no previous run
    # has been performed)   
    if not os.path.exists(uniprot_hhsearch):

        # Inform the user that HHsearch is about to be run
        infostr = \
            f"Running HHsearch for the protein with UniProt ID " \
            f"{uniprot_id}..."
        logger.info(infostr)
        
        # Search profile HHM database(s) with HHsearch
        hhsearch_process = run_process(**hhsearch_params_copy)
    
    # If the output file exists, do not run HHsearch
    else:

        # Inform the user
        warnstr = \
            f"The file {uniprot_hhsearch} was found in " \
            f"{hhsearch_dirname}. SLiMfast will not run HHsearch " \
            f"again."
        logger.warning(warnstr)
        
        # No data will be associated with the HHsearch process,
        # since HHsearch was not run
        hhsearch_process = None

    # Return data about the processes run
    return (hhblits_process, hhsearch_process), uniprot_hhsearch


def run_dssp(pdb_file,
             dssp_params,
             wd):
    """Run DSSP to extract secondary structure information
    from a PDB file.
    """

    # Make a copy of the parameters set to run DSSP
    dssp_params_copy = dict(dssp_params)

    # Path to the directory where DSSP will be run
    dssp_dirname = os.path.join(wd, dssp_params["wd"])

    # Add the directory where DSSP should be run
    dssp_params_copy["wd"] = dssp_dirname

    # Set the path to the DSSP output file
    dssp_output = \
        os.path.join(dssp_dirname,
                     os.path.basename(pdb_file)[:-4] + ".dssp")
                
    # Add the PDB file as input to the options for DSSP
    dssp_params_copy["options"] = \
        {**dssp_params_copy["options"],
         **{"-i" : pdb_file, "-o" : dssp_output}}

    # If the output file does not exist yet (no previous run
    # has been performed)     
    if not os.path.exists(dssp_output):

        # Inform the user that DSSP is about to be run
        infostr = \
            f"Running DSSP on the PDB file {pdb_file}..."
        logger.info(infostr)
        
        # Run DSSP
        dssp_process = run_process(**dssp_params_copy)
    
    # If the output file exists, do not run DSSP
    else:

        # Inform the user
        warnstr = \
            f"The file {dssp_output} was found in " \
            f"{dssp_dirname}. SLiMfast will not run DSSP " \
            f"again."
        logger.warning(warnstr)

        # No data will be associated with the DSSP process,
        # since DSSP was not run
        dssp_process = None

    # Return data about the process run
    return dssp_process, dssp_output