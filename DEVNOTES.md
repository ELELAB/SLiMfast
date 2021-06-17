## Developer notes

`SLiMfast` is divided into three main scripts and some utility scripts.

## Utility scripts

### `defaults.py`

It contains default values used by the pipeline (such as names of columns in the output tables, metadata to be collected for each PDB structure, etc.).

### `disprot_tools.py`

It contains the `DisProtParser` class, implementing a parser for annotations provided in the DisProt database.

DisProt annotations are retrieved through the DisProt API, downloading data from `https://www.disprot.org/api/{uniprot_id}`, where `{uniprot_id}` represents the UniProt ID of the protein whose annotations you are interested in.

### `dssp_tools.py`

It contains the `DSSPParser` class, which implements a parser for DSSP output files.

Only parsing for non-verbose DSSP output files has been implemented so far.

### `hhsearch_tools.py`

It contains the `HHsearchParser` class, implementing a parser for HHsearch output files, and several utility functions to manipulate the data parsed, such as:

* `filter_hits`, which filters the dataframe containing the hits found by HHsearch for homologs of the given protein sequence.
* `sort_hits`, which sorts the dataframe containing the hits.
* `get_aligned_sequence`, which finds for a specific hit the sequence aligned to a specific region of the input UniProt sequence, given the numeric range of such region in the UniProt sequence. This function returns not only the aligned sequence and its starting and ending points in the hit, but also whether the sequence perfectly matches the corresponding one in the input UniProt sequence.

### `pdb_tools.py`

It contains several utility functions that retrieve data from either PDB/mmCIF files or from the PDB/PDBe database, such as:

* `get_uniprot_segments`, which gets which portions of a protein sequence (given its UniProt ID) are covered in a given PDB entry (given its PDB ID). The mapping between the UniProt numbering and the PDB numbering of the given file is retrieved from the PDBe API, at `https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{pdb_id}`, where `{pdb_id}` is the PDB ID of the PDB entry of interest.
* `get_pdb_data`, which gets the content of a PDB file associated to a given PDB ID.
* `get_missing_residues`, which gets residues that are missing in a given PDB entry, given the content of the file associated to the entry, by parsing the REMARK 465 field. See [here](https://www.wwpdb.org/documentation/file-format-content/format33/remarks2.html#REMARK%20465) for more information about the REMARK 465 field. 
* `get_structure_metadata`, which gathers metadata such as experimental technique, resolution, etc. for a given PDB entry, given its PDB ID. Metadata are obtained parsing the mmCIF file associated to the entry (without saving it). Here, mmCIF files are used instead of PDB files since their format is more consistent and less error-prone.

### `run_tools.py`

It contains utility functions to run external tools needed the the pipeline, such as:

* `run_process`, which is a generic function to run an external process.
* `run_hhpred`, which runs the HHpred pipeline to find homologs of a given protein which have a PDB structure available by calling first HHblits and then HHsearch.
* `run_dssp`, which runs DSSP on a given PDB file.

### `uniprot_tools.py`

It includes functions for fetching data from the UniProt database and manipulating them, such as:

* `write_uniprot_fasta`, which fetches the FASTA sequence of a protein, given its UniProt ID.
* `get_pdbs`, which collects the PDB IDs of entries containing a protein of interest, given its UniProt ID.
* `get_pdbs_df`, which creates a dataframe containing information about the PDB structures available for a protein of interest, given the protein's UniProt ID.

### `util.py`

It includes miscellanea utility functions for the pipeline, such as:

* `get_arguments`, which determines which command-line arguments should be passed to all three main scripts and collects their values.

* `read_fasta`, which reads header and sequence from a FASTA file.

* `read_slims_file`, which reads the input file containing the SLiMs predicted to be in a protein sequence.

* `add_structures_metadata`, which takes in input a dataframe containing information about PDB entries and adds metadata about them to the dataframe (collected from the mmCIF files corresponding to those entries).

  Internally, it:

  * Calls `get_structure_metadata` from `pdb_tools.py`.

* `get_structures_data_protein`, which creates a comprehensive dataframe containing information about the PDB structures available for a protein of interest, given the protein's UniProt ID.

  Internally, it:

  * Calls `get_pdbs_df` from `uniprot_tools.py`.
  * Calls `adds_structures_metadata`.

* `get_structures_data_homologs`, which creates a comprehensive dataframe containing information about the PDB structures of the homologs available for a protein of interest, given the protein's UniProt ID.

  Internally, it:

  * Calls `run_hhpred` from `run_tools.py`.
  * Creates a `HHsearchParser` instance from `hhsearch_tools` and call its `parse` method to parse the HHsearch output file.
  * Calls `add_structures_metadata`.
  * Calls  `filter_hits` from `hhsearch_tools` to filter the hits found.
  * Calls `sort_hits` from `hhsearch_tools` to sort the hits found.

* `region_in_protein_segments`, which checks if a specific region of the protein of interest is present in at least one segment of the sequence covered by a given PDB structure.

* `region_in_homolog_segments`, which checks if a specific region of the protein of interest is found in the segments corresponding to such a region in a given homolog of the protein.

  Internally, it:

  * Calls `get_aligned_sequence` from `hhsearch_tools` to find the portion of the homolog corresponding to the region of interest.
  * Calls `get_missing_residues` from `pdb_tools` to find possible missing residues in the PDB structure of the homolog.

* `find_slims`, which look for the SLiMs of interest in the structures of the protein or in those of the homologs.

  Internally, for each PDB structure, it:

  * Calls `get_pdb_data` from `pdb_tools` to fetch structural data for the PDB structure.
  * Calls `region_in_protein_segments` to check if the SLiM is in one of the regions covered by the PDB structure.
  * Calls `region_in_homolog_segments` to check if the SLiM is in one of the regions of the protein covered by the PDB structure of the homolog (called only for PDB structures of homologs, not for those of the protein itself).

* `get_dssp_data`, which runs DSSP and parse its output file.

  Internally, it:

  * Calls `run_dssp` from `run_tools` to run DSSP.
  * Creates a `DSSPParser` instance from `dssp_tools` and calls its `parse` method to parse the DSSP output file.

* `is_region_disordered`, which checks if a SLiM (or, more generally, the region it belongs to) is disordered.

* `get_results`, which checks which SLiMs are in disordered regions and which are not and returns two data frames storing the results for SLiMs found in structures of the protein of interest or of homologs of the protein.

  Internally, for each PDB structure, it:

  * Calls `get_dssp_data` to run DSSP and get information about the secondary structure from the DSSP output file.
  * Calls `is_region_disordered` for each SLiM under investigation to see if it is in a disordered region.

* `get_slims_disprot_annotations`, which collects DisProt annotations about disorder for the SLiMs that could not be found in any structure.

  Internally, it:

  * Creates a `DisProtParser` instance from `disprot_tools` and call its `get_annotations` method to get the annotations from DisProt.

## Main scripts

* `slimfast_get_structures.py`, to get a table containing info about the available PDB structures, given a UniProt ID.
* `slimfast_get_homologs.py`, to get a table containing info about the homologs found for a protein, given its UniProt ID.
* `slimfast_find_slims.py`, to search for SLiMs in the PDB structures contained in a table obtained with either of the previous scripts (or both), and check if they are in a disordered region.