# The UniProt ID of the protein of interest
uniprot_id=$1

# The input CSV file containing the SLiMs to be found
is=slims.csv

# The configuration files for running
crs=run_get_structures.yaml
cra=run_get_alphafold.yaml
crf=run_find_slims.yaml

# The name of the output CSV file containing data about the structures
# of the protein of interest found in the Protein Data Bank
op=${uniprot_id}_structures.csv
oa=${uniprot_id}_alphafold.csv

# The name of the log files
lfs=${uniprot_id}_structures.log
lfa=${uniprot_id}_alphafold.log
lff=${uniprot_id}_slims.log

# Run with verbose logging ('-v'), and log both to the log file
# ('-lf' log_file) and to the standard output ('-ls')
slimfast_get_structures -id $uniprot_id -cr $crs -op $op -lf $lfs -vv -ls

# Run with verbose logging ('-v'), and log both to the log file
# ('-lf' log_file) and to the standard output ('-ls')
slimfast_get_alphafold -id $uniprot_id -cr $cra -oa $oa -lf $lfa -v -ls

# Run with verbose logging ('-v'), and log both to the log file
# ('-lf' log_file) and to the standard output ('-ls')
slimfast_find_slims -id $uniprot_id -cr $crf -is $is -ia $oa -lf $lff -vv -ls

# Run with verbose logging ('-v'), and log both to the log file
# ('-lf' log_file) and to the standard output ('-ls')
slimfast_find_slims -id $uniprot_id -cr $crf -is $is -ip $op -lf $lff -vv -ls
