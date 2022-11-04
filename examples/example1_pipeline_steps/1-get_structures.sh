# UniProt ID of the protein of interest
id=P10415

# Configuration file for running
cr=run_example1.yaml

# Name of the output CSV file containing data about the structures
# of the protein of interest found in the Protein Data Bank
op=${id}_structures.csv

# Name of the log file
lf=${id}_structures.log

# Run with verbose logging (-v), and log both to the log file
# (-lf log_file) and to the standard output (-ls)
slimfast_get_structures -id $id -cr $cr -op $op -lf $lf -v -ls
