# example3_naccess

In this example, the four tools constituting the SLiMfast pipeline (`slimfast_get_structures`, `slimfast_get_alphafold`, `slimfast_get_homologs`, and `slimfast_find_slims`) are run one by one so that each step of the pipeline can be run independently with the NACCESS analysis.

# NACCESS support
If we want to calculate the solvent accessibility with NACCESS, we need to add a column (slim_acc_res) to the csv file with the list of SLiMs to be retrieved. The residues for which we want to calculate the SASA need to be ':' separated (i.e., 1280:1283)

#RUN
First, we want to get the IDs of the protein structures of the protein of interest available in the Protein Data Bank and some associated metadata. We use the `1-get_structures.sh` script, which conveniently wraps the call to `slimfast_get_structures`.

```shell
sh 1-get_structures.sh
```
Then, we want to retrieve the IDs of the protein structures of the protein of interest available in the Alphafold Protein Structure Database and associated metadata. Here, we run the `2-get_alphafold.sh` script, which runs `slimfast_get_alphafold` on our behalf.

```shell
sh 2-get_alphafold.sh
```

Furthermore, we want to see if the protein of interest has some homologs whose structures are available in the Protein Data Bank if no structures of the protein itself are found either in the Protein Data Bank or the Alphafold Protein Structure Database. To do this, we use the `3-get_homologs.sh` script, devised to run the `slimfast_get_homologs` executable.

```shell
sh 3-get_homologs,sh
```

Finally, we want to see if we find our SLiMs of interest in either the structures of the protein of interest or in those of its homologs. Therefore, we take advantage of the `4-find_slims.sh` script, which calls `slimfast_find_slims`.

```shell
sh 4-find_slims.sh
```
