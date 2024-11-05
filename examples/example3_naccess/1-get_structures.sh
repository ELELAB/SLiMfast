# The UniProt ID of the protein of interest
id=Q9BQS8

# The configuration file for running
cr=run_get_structures.yaml

# The name of the output CSV file containing data about the structures
# of the protein of interest found in the Protein Data Bank
op=${id}_structures.csv

# The name of the log file
lf=${id}_structures.log

# Run with verbose logging ('-v'), and log both to the log file
# ('-lf' log_file) and to the standard output ('-ls')
slimfast_get_structures -id $id -cr $cr -op $op -lf $lf -vv -ls
