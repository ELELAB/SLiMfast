# In this example, the four tools constituting the SLiMfast pipeline
# (slimfast_get_structures, slimfast_get_alphafold,
# slimfast_get_homologs, and slimfast_find_slims) are run one by one,
# so that each step of the pipeline can be run on its own.


# First, we want to get the IDs, chain IDs, and some associated
# metadata of the structures of the protein of interest available
# in the Protein Data Bank.

# For this, we use the 1-get_structures.sh script, which conveniently
# wraps the call to slimfast_get_structures.

sh 1-get_structures.sh


# Then, we want to retrieve the IDs and some associated metadata
# of the structures of the protein of interest available in the
# AlphaFold Protein Structure Database.

# To achieve this goal, we run the 2-get_alphafold.sh script, which
# runs the slimfast_get_alphafold executable

sh 2-get_alphafold.sh


# Then, we also want to see if the protein of interest has some
# homologs whose structures are available in the Protein Data Bank,
# in case no structures of the protein are found in either the
# Protein Data Bank or the AlphaFold Protein Structure Database.

# To do this, we use the 3-get_homologs.sh script, devised to run
# the slimfast_get_homologs tool.

sh 3-get_homologs.sh


# Finally, we want to see if we find our SLiMs of interest in either
# the structures of the protein of interest found in the Protein Data
# Bank/AlphaFold Protein Structure Database or in the structures
# of the homologs found in the Protein Data Bank.

# Therefore, we take advantage of the 4-get_homologs.sh script, which
# calls the slimfast_find_slims tool.

sh 4-find_slims.sh

