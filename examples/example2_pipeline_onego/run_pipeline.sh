# File containing the UniProt IDs of the proteins of interest
# and the SLiMs to be found in the protein sequences
ids=ids_and_slims.csv

# Configuration file for running
cr=run_example2.yaml

# Number of processes to be used when running
n=8

# Run with verbose logging (-v), and log both to the log file
# (-lf log_file) and to the standard output (-ls). Run all steps
# of the SLiMfast pipeline in one go (-dgs, -dga, -dgh, -dfs).
slimfast_run_pipeline -dgs -dga -dgh -dfs -ids $ids -cr $cr -ls -v -n $n
