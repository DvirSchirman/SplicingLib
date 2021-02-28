# Introduction
This repository contanis the data analysis and alignment code for the paper doi.org/10.1101/2020.04.20.050609

# Data analysis
all the data analysis code are found under analysis/Source/.
All scripts are written in Matlab (R2019b), except for the gradient boosting model analysis (Figure 5 + Figure 6E) that are written in python (3.8).
The scripts are organized according to figure numbers in the paper.
Scripts that are not named by figure number are used for processing of data files used by other scripts.

# Mapping NGS data
The scripts under "Mapping NGS data" are used to map the raw sequencing reads to the library database and to compute from it splicing efficiency values.
The scripts are written to run on a LSF based cluster, and use bioinformatic tools installed on the cluster.
All the scripts are written in Matlab (R2019b). The scripts are numbered in the order they should be run, and unnumbered are functions that are called by another script.
The raw sequencing files can be found in the Sequence Read Archive (SRA), Accession number PRJNA631112.
For running the maping pipeline, download SplicingLib1 fastq files from SRA and put the files in "mapping NGS data/Data/Samples" folder.

## Tools required:
* PEAR
* Cutadapt
* vsearch


