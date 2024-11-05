# The UniProt ID of the protein of interest
id=Q9BQS8

# The configuration file for running
cr=run_get_alphafold.yaml

# The name of the output CSV file containing data about the structures
# of the protein of interest found in the AlphaFold Protein Structure
# Database
oa=${id}_alphafold.csv

# The name of the log file
lf=${id}_alphafold.log

# Run with verbose logging ('-v'), and log both to the log file
# ('-lf' log_file) and to the standard output ('-ls')
slimfast_get_alphafold -id $id -cr $cr -oa $oa -lf $lf -v -ls
