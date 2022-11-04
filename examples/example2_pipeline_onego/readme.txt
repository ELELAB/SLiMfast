# In this example, the four tools constituting the SLiMfast pipeline
# (slimfast_get_structures, slimfast_get_alphafold,
# slimfast_get_homologs, and slimfast_find_slims) are run in one go
# for a set of proteins and their SLiMs of interest, thanks to the
# slimfast_run_pipeline tool.

# To run the pipeline, it is sufficient to launch the run_pipeline.sh
# script where input and output files are defined together with the
# steps of the SLiMfast pipeline to run (in this example, all of them).

sh run_pipeline.sh
