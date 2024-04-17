# The UniProt ID of the protein of interest
id=P10415

# The configuration file for running
cr=run_get_homologs.yaml

# The name of the output CSV file containing data about the structures
# of the homologs of the protein of interest found in the Protein Data
# Bank
oh=${id}_homologs.csv

# The name of the output JSON file containing the alignments of the
# portions of the sequence of each homolog with the corresponding
# portions of the sequence of the protein of interest
oal=${id}_homologs.json

# The name of the log file
lf=${id}_homologs.log

# Run with verbose logging ('-v'), and log both to the log file
# ('-lf' log_file) and to the standard output ('-ls')
slimfast_get_homologs -id $id -cr $cr -oh $oh -oal $oal -lf $lf -v -ls
