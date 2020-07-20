#!/usr/bin/env bash
##############################################################################
# Author: Joe Colgan                   Program: extract_greatest_value_per_gene.sh
#
# Date: 28/02/2019
#
# Purpose:
# The script takes a BED file as input.
# The first step involves intersecting the BED file that contains individual
# gene coordinates generated from a ensembl gff3 file (release 42) for the 
# buff-tailed bumblebee with a second BED file containing SNP positions and an 
# assigned standardised iHS or nSL value.
# The intersection is sorted by highest iHS or nSL value and output.
# The output file is a tab-delimited text file containing two columns:
# 1) Locus ID
# 2) iHS or nSL value for that gene
#
##############################################################################

## Module load programs required:
module load bedtools

## Take input as arguments from the command line:
input=$1 ## The first argument on the command line will be taken as input
output=$2

## Check input arguments exist. 
## If not, print usage.
if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "Usage: extract_greatest_value_per_gene.sh input.bed output.txt"
fi

## Extract chromosome name:
chrom_name="$(echo "$input" | cut -d '-' -f 1 )"
echo "$chrom_name"

## Extract locus name:
locus_name="$(echo "$input" | cut -d '-' -f 4 | cut -d '.' -f 1 )"
echo "$locus_name"

## Print to console:
echo "Running" 

## Intersect input file with a reference BED file containing iHS or nSL values:
value="$(intersectBed -a ../"$chrom_name".ihs.out.100bins.norm.abs.bed.tmp -b "$input" | \
            sort -k4,4gr | \
            head -n 1 | \
            awk '{ print $4 }')"

## Print to output:
echo "$locus_name" "$value" | sed 's/ /\t/g' >> "$output"

## Print to console:
echo "Complete"
