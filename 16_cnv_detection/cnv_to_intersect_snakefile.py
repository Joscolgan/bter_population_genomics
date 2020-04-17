#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: CNV_to_intersect.py
#
# Date: 10/01/2016
#
##############################################################################
# Import modules
import os.path
import time
import sys
import subprocess

from helper_functions import *

# Script information
# This script takes one tab-delimited file (output by the software program, CNVnator)
# as input, quality filters and reformats into BED format. The resulting BED file is
# compared against a BED file containg gene annotations for species of interest to
# identify intersection between high quality duplication events and potential coding
# sequences. The final output is a BED file containing a matrix populated by:
# Columns: Duplication
# Samples: Sample name
#
# To achieve final output, the Snakefile contains custom-definedr rules (see below) that
# outline commands to execute sequentially to take custom-defined input(s) and generate
# final output (as defined in rule all).
#
# For this specific script, rules are defined as follows:
# rule all:
# - Defines the expected output of the Snakefile.
# rule filter_by_pvalue:
# - For each input, CNV calls will be filtered by p-value.
# rule filter_by_mapping:
# - For each input, CNV calls will be filtered by mapping quality.
# rule filter_by_size:
# - For each input, CNV calls will be filtered by size.
# rule convert_to_BED:
# - For each sample, tab-delimited input containing filtered CNV calls
# will be reformatted into BED file format using a custom-designed script
# 'cnv_to_bed.py'.
# rule extract_duplications:
# - For each input, duplication annotated CNVs will be paresed.
# rule intersect_beds:
# - For each input, intersects are performed for co-ordinates between each
# sample and BED file containing gene information for species of interest and output
# as a BED file containing co-ordinates contained within both sample and reference.

##############################################################################
# Sample information
##############################################################################
# Bumblebee (Bombus terrestris) males were collected summer 2014
# Sample sites (n=26) were from across the UK
# Each individual and site were assigned unique identifiers
# Example: '2014_Bter_P_D_14_260_head'
# Explanation:{year_collected}_{species}_{site_type}_{sex}_{site_number}_{tube_number}_{tissue_type}
# species: Bter = Bombus terrestris
# site_type: P = Pastoral, M = Mixed, A = Arable
# sex: D = Male (drone)

##############################################################################
#
# Prior to use
#
##############################################################################

# To run read_filtering_snakefile.py:
#  1. Download and install the following software:
#  bedtools2 (https://github.com/arq5x/bedtools2.git)
#  OPTIONAL: BEDOPS (http://bedops.readthedocs.org/en/latest/) for the conversion of reference gff
#  file into BED format

#  2. Assign global variables for use within specific rules
#   Please see section below for further information on variable to be assigned

#  3. Assignment of wildcards to be used within the rules

#  4. Define all the input and outputs of each rule with respect to the above defined wildcards

#  5. Each rule will take assigned input, global variables and generate defined outputs.

#  6. Input data should be formatted in the context of defined wildcards: {sample}.{pair}.fastq
#       For example: 2014_Bter_P_D_14_260_head.R1.fastq

#  7. Make a text.file called 'sample_list.txt' and put in same directory as Snakfile.\
#    Populate 'sample_list.txt' with names of samples to be analysed.
#       For example:
#       2014_Bter_P_D_14_260_head
#       2014_Blap_P_D_21_412_thorax
#       2014_Bpas_A_D_20_391_thorax

##############################################################################
#
# Running Snakemake
#
##############################################################################

#  To run Snakefile present within current directory
#  snakemake -s <Snakefile_name>.py -p -j <number of cores>
#  '-s' option is required if Snakefile has customs named
#  '-p' prints shell commands to console during execution
#  '-dag | dot -Tpdf > dag.pdf' gives a nice pretty graph of what will happen

##############################################################################
# Assign global variables for use in rules (see below)
##############################################################################

# For rule filter_by_pvalue, assign value for p-value cut-off.
P_VALUE = 0.05

# For rule filter_by_mapping, assign value for mapping quality cut-off.
MAPQUAL = 0.5

# For rule filter_by_size, assign value for CNV size cut-off.
SIZE = 1000

# For rule intersect_genes, assign path for BED file containing gene information
# extracted from Bombus reference GFF file.
REF_BED  = "../../../../2015-12-24_subsample_BAM_test/GCF_000214255.1_Bter_1.0_genomic.gff.gene.bed"

##############################################################################
# Assignment of wildcards to be used within rules
##############################################################################

# Open file and read in contents - one sample per line

SAMPLES = open('samples_list.txt').read().splitlines()
for sample in SAMPLES:
    if sample.find('.') != -1:
        raise IOError("Names should use '_', not '.'. Problem with: %s" % sample)

##############################################################################
# Specify all input/output files in terms of sample wildcards
##############################################################################

# Assign path for subsample data output
RAW_DATA             = "../{samples}.txt"

# Assign path for p-value filtered data output
PVALUE_FILTERED_DATA = "temp/01_filtered_data/{samples}.pvalue_filtered.txt"

# Assign path for mapping quality filtered data output
MAPQ_FILTERED_DATA   = "temp/01_filtered_data/{samples}.mapq_filtered.txt"

# Assign path for size filtered data output
SIZE_FILTERED_DATA   = "temp/01_filtered_data/{samples}.size_filtered.txt"

# Assign path for BED formatted data output
BED_DATA             = "temp/02_BED_converted/{samples}.bed"

# Assign path for BED file containing duplication information output
DUPLICATED_DATA      = "temp/03_duplicated_data/{samples}.duplicated.bed"

# Assign path for intersected data output
INTERSECT_DATA       = "temp/04_intersected/{samples}.intersect_with_GFF.bed"

#############################################################################
# Define binaries in context of path relative to Snakefile
##############################################################################
# binaries
# Align lines of code using 'Assign Align' using cmd+shift+p and selecting 'align'
# Create dictionaries for directories and  tools
dirs  = {}
dirs['project'] = os.path.abspath('../../../../../../')
dirs['src']     = os.path.join(dirs['project'], 'src')

tools = {}
tools['intersect'] = os.path.join(dirs['src'], 'bedtools2/bin/intersectBed')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################

# First rule is list the final output
rule all:
    input: expand(INTERSECT_DATA,\
                  samples = SAMPLES)

# Filter CNV calls by p-value (p<0.05) and output
rule filter_by_pvalue:
    input: RAW_DATA
    output: PVALUE_FILTERED_DATA
    run:
        check_files_arent_empty(input)
        shell("awk '$5 < {P_VALUE} {{print $0}}' {input} > {output}")
        check_files_arent_empty(output)

# Filter CNV calls by mapping quality (<0.50) and output
rule filter_by_mapping:
    input: PVALUE_FILTERED_DATA
    output: MAPQ_FILTERED_DATA
    run:
        check_files_arent_empty(input)
        shell("awk '$9 < {MAPQUAL} {{print $0}}' {input} > {output}")
        check_files_arent_empty(output)

# Filter CNV calls by size (>1000 bases) and output
rule filter_by_size:
    input: MAPQ_FILTERED_DATA
    output: SIZE_FILTERED_DATA
    run:
        check_files_arent_empty(input)
        shell("awk '$3 > {SIZE} {{print $0}}' {input} > {output}")
        check_files_arent_empty(output)

# Convert CNV calls to BED format and output
rule convert_to_BED:
    input: SIZE_FILTERED_DATA
    output: BED_DATA
    run:
        check_files_arent_empty(input)
        shell("python cnv_to_bed.py {input} {output}")
        check_files_arent_empty(output)

# Subset CNV calls annotated as 'duplication' and output
rule extract_duplications:
    input: BED_DATA
    output: DUPLICATED_DATA
    run:
        check_files_arent_empty(input)
        shell("grep 'duplication' {input} > {output}")
        check_files_arent_empty(output)

# Intersect CNV BED files with GFF BED file and output
rule intersect_genes:
    input: DUPLICATED_DATA
    output: INTERSECT_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[intersect]} -a {REF_BED} -b {input} > {output}")
        check_files_arent_empty(output)
