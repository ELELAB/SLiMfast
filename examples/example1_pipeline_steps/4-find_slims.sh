# The UniProt ID of the protein of interest
id=P10415

# The configuration file for running
cr=run_find_slims.yaml

# The input CSV file containing the SLiMs to be found
is=slims.csv

# The input CSV file containing data about the structures of the
# protein of interest found in the Protein Data Bank
ip=${id}_structures.csv

# The input CSV file containing data about the structures of the
# protein of interest found in the AlphaFold Protein Structure
# Database
ia=${id}_alphafold.csv

# The input CSV file containing data about the structures of the
# homologs of the protein of interest found in the Protein Data
# Bank
ih=${id}_homologs.csv

# The input JSON file containing the alignments of the portions of
# the sequence of each homolog with the corresponding portions of
# the sequence of the protein of interest
ial=${id}_homologs.json

# The name of the log file
lf=${id}_slims.log

# Run with verbose logging ('-v'), and log both to the log file
# ('-lf' log_file) and to the standard output ('-ls')
slimfast_find_slims -id $id -cr $cr -is $is -ip $ip -ia $ia -ih $ih -ial $ial -lf $lf -v -ls
