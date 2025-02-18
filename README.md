# Introduction
*PSTAR-MSA* (Pairwise STructural Alignment Reference for Multiple Sequence Alignment) is an MSA benchmarking tool that aims to generalize the comparative scoring of multiple sequence alignment (MSA) algorithms. This is tried to be achieved by using pairwise structural alignments as reference data for each pair of aligned sequences in a given MSA. Tertiary protein structure data is sourced from the *"Protein Data Bank"* (PDB) and is planned to be expanded by high-certainty protein structure predictors like *AlphaFold*. Novel measure indices like the **PSTAR-score**, which is the cumulative pairwise structure alignment equivalent of the **sum-of-pairs-score** (SPS), have been itroduced in order to evaluate the similarity between MSA algorithms.
Usually MSA algorithm quality is measured through the use of manually curated reference alignment databases such as *HOMSTRAD*, or *HomFam* for that matter. However, all of these databases suffer from limited data availability and it is not clear, if their status as references is legitimate. Because of this limited data availability it is likely that algorithms perform worse in general than they perform on the small select data set.
An in-depth explanation of the scoring and some first results can be found [<ins>here</ins>](link to pdf). (Link not available yet.)

# Installation and setup
**Linux** environment required.
List of dependencies:
```
sudo apt-get install rsync
sudo apt-get install python3.6
sudo apt install python3-pip
pip install argparse DateTime pandas requests rcsbsearchapi biopython matplotlib
```
install [DaliLite.v5](http://ekhidna2.biocenter.helsinki.fi/dali/) from http://ekhidna2.biocenter.helsinki.fi/dali/
install [DIAMOND](https://github.com/bbuchfink/diamond) from https://github.com/bbuchfink/diamond
download **PSTAR-MSA.py** from the master branch.

In case you're not in a hurry, you are done now. Downloading the PDB database and setting up the DIAMOND database can take up to a few hours. You don't need to create the specified directories manually, but you can. You can start the process via
```
python3 -m PSTAR-MSA --syncdb --diamondfile <path_to_DIAMOND_exe_file> --locpdbdb <dir_path__to_create_local_pdb_copy_in> --outputdir <dir_path_to_place_DIAMOND_output_db_files_in> --verbose
```
### **Waiting is recommended to ensure flawless database generation.**

However, you can download [these files](link to zip), manually create a directory you can later refer to under **--outputdir** of the **--syncdb** option, then extract the two files in the zip archive into this directory and launch the above command with these two files already in place and specifying the directory name under the **--outputdir** flag. The generation of the FASTA file based on the PDB copy is then skipped. But it can happen, if entries were removed from the PDB that manual deletions in the FASTA file and manual subsequent DIAMOND database generation have to be carried out.

# 
