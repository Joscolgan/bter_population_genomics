#!/bin/env bash
##################################################################################################
##
## Author: Joe Colgan                                                       name: run_kallisto.sh
##
## Date:
## 
##
## Purpose:
## This script the following commands:
## - Index of input cDNA file for pseudoalignment.
## - Pseudoaligmment of individual single-end FASTQ files against
##   kallisto-indexed cDNA sequences.  
##
## The script takes compressed (gzip) FASTQ files (one per sample) and results in the generation of 
## three files per sample. 
## 1) abundances.h5: HDF5 binary file containing run info, abundance estimates, bootstrap estimates,
##    and transcript length information.
## 2) abundances.tsv: Plaintext format of abundance estimates. Does not contain boostrap estimates. 
## 3) run_info.json: A json file containing information on the run.  
##
##################################################################################################

## Take input arguments from the command line:
reference=$1
input=$2
output=$3

if [[ $# -eq 0 ]] ; then
    echo 'No arguments provided'
    echo 'Usage: ./run_kallisto.sh index_name path_to_input_sequences output_directory'
    exit 1
fi

## Run kallisto:
kallisto quant \
-i "$reference"  \
--output-dir="$output" \
--threads=40 \
-b 100 \
--bias \
--single \
--fragment-length=300 --sd=20 \
"$input"_pass.fastq.gz
