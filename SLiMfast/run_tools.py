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

    # Create an empty list of arguments
    args = []

    # For each option, value pair
    for opt, val in options.items():

        # If the value is True
        if val is True:

            # Simply append the option (it is a flag without value)
            args.append(str(opt))

        # If the value is None
        elif val == None:

            # Skip the pair altogether
            continue

        # Otherwise
        else:

            # Append both the option and its value to the list
            # of arguments
            args.extend([str(opt), str(val)])

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
            "pid" : popen.pid,
            "returncode" : popen.returncode}


def run_pipeline_step(step,
                      wd,
                      uniprot_id,
                      input_files,
                      output_files,
                      configrun_file,
                      log_file,
                      v,
                      vv):
    """Run a step of the SLiMfast pipeline.
    """

    # Get the name of the log file
    log_file = f"{uniprot_id}_{log_file}"

    # Set the common options
    options = \
        {"--uniprot-id" : uniprot_id,
         "--configrun-file" : configrun_file,
         "--log-file" : log_file}

    # If verbose logging mode has been requested 
    if v:

        # Update the options
        options.update({"-v" : True})

    # If debug mode has been requested
    if vv:

        # Update the options
        options.update({"-vv" : True})

    # If the step consists in running slimfast_get_structures
    if step == "get_structures":

        # Set the executable
        executable = "slimfast_get_structures"

        # Update the options
        options.update(\
            {"--output-protein-file" : output_files[0]})

    # If the step consists in running slimfast_get_alphafold
    elif step == "get_alphafold":

        # Set the executable
        executable = "slimfast_get_alphafold"

        # Update the options
        options.update(\
            {"--output-alphafold-file" : output_files[0]})

    # If the step consists in running slimfast_get_homologs
    elif step == "get_homologs":

        # Set the executable
        executable = "slimfast_get_homologs"

        # Update the options
        options.update(\
            {"--output-homologs-file" : output_files[0],
             "--output-alignments-file" : output_files[1]})

    # If the step consists in running slimfast_find_slims
    elif step == "find_slims":

        # Get the executable
        executable = "slimfast_find_slims"

        # Unpack the input files
        slims_input_file, protein_file, alphafold_file, \
            homologs_file, alignments_file = input_files

        # Unpack the output files
        slims_protein_file, slims_alphafold_file, \
            slims_homologs_file, slims_disprot_file = output_files

        # Add the SLiMs' input file
        options.update({"--input-slims-file" : slims_input_file})

        # If the file containing the structures of the protein
        # of interest found in the Protein Data Bank is not None
        if protein_file is not None:

            # Get the full path of the file
            protein_file = os.path.join(wd, protein_file)

            # If the file exists
            if os.path.exists(protein_file):

                # Add the corresponding options
                options.update(\
                    {"--input-protein-file" : \
                        protein_file,
                     "--output-slims-protein-file" : \
                        slims_protein_file})

        # If the file containing the structures of the protein
        # of interest found in the AlphaFold Protein Structure
        # Database is not None
        if alphafold_file is not None:

            # Get the full path of the file
            alphafold_file = os.path.join(wd, alphafold_file)

            # If the file exists
            if os.path.exists(alphafold_file):

                # Add the corresponding options
                options.update(\
                    {"--input-alphafold-file" : \
                        alphafold_file,
                     "--output-slims-alphafold-file" : \
                        slims_alphafold_file})

        # If the file containing the structures of the homologs
        # of the protein of interest found in the Protein Data
        # Bank and the file containing the alignments
        # are not None
        if homologs_file is not None and alignments_file is not None:

            # Get the full path of the files
            homologs_file = os.path.join(wd, homologs_file)
            alignments_file = os.path.join(wd, alignments_file)

            # If the files exists
            if os.path.exists(homologs_file) \
            and os.path.exists(alignments_file):

                # Add the corresponding options
                options.update(\
                    {"--input-homologs-file" : \
                        homologs_file,
                     "--input-alignments-file" : \
                        alignments_file,
                     "--output-slims-homologs-file" : \
                        slims_homologs_file})

    # Launch the process
    process = \
        run_process(executable = executable,
                    options = options, 
                    wd = wd)

    # Return the output files
    return output_files


def update_hhblits_params(hhblits_params,
                          wd,
                          uniprot_fasta,
                          uniprot_id):
    """Update the HHblits parameters before running.
    """

    # Optional output files that can be generated by HHblits
    OPT_OUT_FILES = \
        (
         # Write results in standard format to file
         "-o",
         # Write result MSA of significant matches in PSI-BLAST format
         "-opsi",
         # Write HHM file for result MSA of significant matches
         "-ohhm", 
         # Write MSAs in A3M format after each iteration
         "-oalis",
         # Write result in tabular BLAST format (compatible
         # with -m 8 or -outfmt 6 output)
         "-blasttab",
         # Write pairwise alignments in FASTA format
         "-Ofas",
         # Write pairwise alignments in A2M format
         "-Oa2m",
         # Write pairwise alignments in A3M format
         "-Oa3m",
         # Write scores for all pairwise comparisons to file 
         "-scores",
         # Write all alignments in tabular layout to file
         "-atab",
         )

    # Make a copy of the parameters set to run HHblits
    # (so that the original dictionary is not modified)
    hhblits_params_copy = dict(hhblits_params)

    # Set the path to the directory where HHblits will be run
    hhblits_dirname = os.path.join(wd, hhblits_params["wd"])

    # Set the path to the directory where HHblits will be run
    hhblits_params_copy["wd"] = hhblits_dirname

    # For each of the options provided to run HHblits
    for opt, val in hhblits_params["options"].items():

        # If the option is one of the possible optional
        # output files
        if opt in OPT_OUT_FILES:

            # Set the path for the output file
            new_val = os.path.join(hhblits_dirname, val)

            # Add it to the updated dictionary of parameters
            hhblits_params_copy["options"].update({opt : new_val})

    # If the option -o has not been passed
    if "-o" not in list(hhblits_params_copy["options"].keys()):

        # Pass it and give the correct path to the file,
        # otherwise the file willbe written anyway, but in
        # the wrong location
        hhblits_params_copy["options"].update(\
            {"-o" : os.path.join(hhblits_dirname,
                                 uniprot_id + "_hhblits.hhr")})

    # Set the path to the output A3M file
    uniprot_a3m = \
        os.path.join(hhblits_dirname, uniprot_id + ".a3m")

    # Set the path to the log file
    hhblits_log = \
        os.path.join(hhblits_dirname, hhblits_params["outfile"])

    # Set the path to the log file
    hhblits_params_copy["outfile"] = hhblits_log

    # Add the FASTA file as input and the A3M file as output
    # to the options of HHblits
    hhblits_params_copy["options"].update(\
        {"-i" : uniprot_fasta, "-oa3m" : uniprot_a3m})

    # Return the updated dictionary of parameters
    return hhblits_params_copy


def update_hhsearch_params(hhsearch_params,
                           wd,
                           uniprot_a3m,
                           uniprot_id):

    # Optional output files that can be generated by HHsearch
    OPT_OUT_FILES = \
        (
         # Write result MSA with significant matches in A3M format
         "-oa3m",
         # Write result MSA of significant matches in PSI-BLAST format
         "-opsi",
         # Write HHM file for result MSA of significant matches
         "-ohhm", 
         # Write result in tabular BLAST format (compatible
         # with -m 8 or -outfmt 6 output)
         "-blasttab",
         # Write pairwise alignments in FASTA format
         "-Ofas",
         # Write pairwise alignments in A2M format
         "-Oa2m",
         # Write pairwise alignments in A3M format
         "-Oa3m",
         # Write scores for all pairwise comparisons to file 
         "-scores",
         # Write all alignments in tabular layout to file
         "-atab",
         )
    
    # Make a copy of the parameters set to run HHsearch
    hhsearch_params_copy = dict(hhsearch_params)

    # Set the path to the directory where HHsearch will be run
    hhsearch_dirname = os.path.join(wd, hhsearch_params["wd"])

    # Set the path to the directory where HHsearch will be run
    hhsearch_params_copy["wd"] = hhsearch_dirname

    # For each of the options provided to run HHblits
    for opt, val in hhsearch_params["options"].items():

        # If the option is one of the possible optional
        # output files
        if opt in OPT_OUT_FILES:

            # Set the path for the output file
            new_val = os.path.join(hhsearch_dirname, val)

            # Add it to the updated dictionary of parameters
            hhsearch_params_copy["options"].update({opt : new_val})

    # Path to the HHsearch output
    uniprot_hhsearch = \
        os.path.join(hhsearch_dirname, uniprot_id + ".hhr")

    # Set the path to the log file
    hhsearch_log = \
        os.path.join(hhsearch_dirname, hhsearch_params["outfile"])

    # Set the path to the log file
    hhsearch_params_copy["outfile"] = hhsearch_log

    # Add the A3M file as input and the output
    hhsearch_params_copy["options"].update(\
        {"-i" : uniprot_a3m, "-o" : uniprot_hhsearch})

    # Return the updated dictionary of parameters
    return hhsearch_params_copy


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


    # Get an updated dictionary of parameters to run HHblits
    hhblits_params_updated = \
        update_hhblits_params(hhblits_params = hhblits_params,
                              wd = wd,
                              uniprot_fasta = uniprot_fasta,
                              uniprot_id = uniprot_id)

    # Get the name of the output file
    uniprot_a3m = hhblits_params_updated["options"]["-oa3m"]

    # If the output file does not exist yet (no previous run
    # has been performed)
    if not os.path.exists(uniprot_a3m):

        # Inform the user that HHblits is about to be run
        infostr = \
            f"Running HHblits for the protein with UniProt ID " \
            f"{uniprot_id}..."
        logger.info(infostr)
        
        # Run HHblits for the generation of the query MSA
        hhblits_process = run_process(**hhblits_params_updated)
    
    # If the output file exists, do not run HHblits
    else:

        # Inform the user
        warnstr = \
            f"The file {uniprot_a3m} was found for the protein with " \
            f"UniProt ID {uniprot_id}. SLiMfast will use it as " \
            f"input to HHsearch, without running HHblits again."
        logger.warning(warnstr)

        # No data will be associated with the HHblits process in this
        # case, since HHblits was not run
        hhblits_process = None


    #-------------------------- Run HHsearch -------------------------#


    # Get an updated dictionary of parameters to run HHsearch
    hhsearch_params_updated = \
        update_hhsearch_params(hhsearch_params = hhsearch_params,
                               wd = wd,
                               uniprot_a3m = uniprot_a3m,
                               uniprot_id = uniprot_id)

    # Get the name of the output file
    uniprot_hhsearch = hhsearch_params_updated["options"]["-o"]

    # If the output file does not exist yet (no previous run
    # has been performed)   
    if not os.path.exists(uniprot_hhsearch):

        # Inform the user that HHsearch is about to be run
        infostr = \
            f"Running HHsearch for the protein with UniProt ID " \
            f"{uniprot_id}..."
        logger.info(infostr)
        
        # Search profile HHM database(s) with HHsearch
        hhsearch_process = run_process(**hhsearch_params_updated)
    
    # If the output file exists, do not run HHsearch
    else:

        # Inform the user
        warnstr = \
            f"The file {uniprot_hhsearch} was found for the " \
            f"protein with UniProt ID {uniprot_id}. SLiMfast will " \
            f"not run HHsearch again."
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
