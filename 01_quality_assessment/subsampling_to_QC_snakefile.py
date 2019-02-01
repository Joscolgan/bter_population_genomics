#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: subsampling_to_QC_snakefile.py
#
# Date: 05/03/2015
#
##############################################################################
# Import modules
import os.path
import time
import sys

from helper_functions import *

# This script takes two fastq files (i.e. pairs) per sample as input and performs
# a quality assessment of a subsampled dataset for one or more samples.
# Final output of script is a QC report in two formats (pdf and html).
#
# To achieve final output, the Snakefile contains user-defined rules (see below) that
# outline commands to execute sequentially to take user-defined input and generate
# final output (as defined by user in rule all).
#
# For this script, three rules are defined.
# 1) rule subsample_each_pair:
#    - For each input sample, a specific number of random sequences are subsampled from each pair.
# 2) rule combine_each_pair:
#  - For each input sample, subsampled sequences are combined per pair.
# 3) rule fastqc_each_pair:
#  - Using the combined subsampled sequences as input, quality assessment is performed\
#   for each pair. This step results in the required output defined in rule all.

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
#   seqtk (https://github.com/lh3/seqtk)
#   fastqc (https://github.com/dib-lab/khmer)
#
#  2. Ensure helper_functions.py is within the same directory of the Snakefile.
#
#  3. Assign global variables for use within specific rules
#   Please see section below for further information on variable to be assigned
#
#  4. Assignment of wildcards to be used within the rules
#
#  5. Define all the input and outputs of each rule with respect to the above defined wildcards
#
#  6. Each rule will take assigned input, global variables and generate defined outputs.
#
#  7. Input data should be formatted in the context of defined wildcards: {sample}.{pair}.fastq
#       For example: 2014_Bter_P_D_14_260_head.R1.fastq
#
#  8. Make a text.file called 'sample_list.txt' and put in same directory as Snakfile.\
#    Populate 'sample_list.txt' with names of samples to be analysed.
#       For example:
#       2014_Bter_P_D_14_260_head
#       2014_Blap_P_D_21_412_thorax
#       2014_Bpas_A_D_20_391_thorax

##############################################################################
# Running Snakemake
##############################################################################

# To run Snakefile present within current directory
# snakemake -s <Snakefile_name>.py -p -j <number of cores>
# '-s' option is required if Snakefile is not named Snakefile (e.g. subsampling_to_QC_snakefile.py)
# '-p' prints shell commands to console during execution of each rule.
# '-j' number of cores. Snakemake can run jobs in parallel.

##############################################################################
# Assign global variables for use in rules (see below)
##############################################################################

# For rule subsample_each_pair, assign path for working data
RAW_DATA  = "data/{samples}.{pair}.fastq"

# For rule subsample_each_pair, specify a number (INT) of sequences to subsample
SUBSAMPLE = 100000

# For rule fastqc_each_pair; Snakemake will use min(MAX_THREADS, --cores)
MAX_THREADS = 40

##############################################################################
# Assignment of wildcards to be used within rules
##############################################################################

# Open file and read in contents - one sample per line

SAMPLES = open('sample_list.txt').read().splitlines()
for sample in SAMPLES:
    if sample.find('.') != -1:
        raise IOError("Names should use '_', not '.'. Problem with: %s" % sample)

# Assign pair information
PAIR = ["R1", "R2"]

##############################################################################
# Specify all input/output files in terms of sample wildcards
##############################################################################

# Assign path for subsample data output
RANDOM_SUBSAMPLE_DATA = "temp/01_subsample/{samples}.{pair}.{subsample}.fastq"

# Combine subsampled data and output
COMBINED_DATA         = "temp/02_combined/combined.{pair}.{subsample}.fastq"

# Assign directory to contain information from fastqc analysis
FASTQC_DATA           = "temp/02_combined/combined.{pair}.{subsample}_fastqc.html",\
                        "temp/02_combined/combined.{pair}.{subsample}_fastqc.zip"

#############################################################################
# Define binaries in context of path relative to Snakefile
##############################################################################
# binaries
# Align lines of code using 'Assign Align' using cmd+shift+p and selecting 'align'
# Create dictionaries for directories and  tools
dirs  = {}
dirs['project'] = os.path.abspath('../../')
dirs['src']     = os.path.join(dirs['project'], 'src')

tools = {}
tools['seqtk']  = os.path.join(dirs['src'], 'seqtk/seqtk')
tools['fastqc'] = os.path.join(dirs['src'], 'FastQC/fastqc')
check_tools(tools)

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################

# First rule is list the final output
rule all:
    input: expand(FASTQC_DATA, pair=PAIR, subsample=SUBSAMPLE)

# Subsample sequences per pair per sample
rule subsample_each_pair:
    input: expand("data/{{samples}}.{pair}.fastq", pair=PAIR)
    output: expand("temp/01_subsample/{{samples}}.{pair}.{subsample}.fastq",
                   pair=PAIR,
                   subsample=SUBSAMPLE)
    run:
        check_identical_fastq_headers(input[0], input[1])
        check_files_arent_empty(input)
        seed = int(time.time() * 10000)
        shell("{tools[seqtk]} sample -s {seed} {input[0]} {SUBSAMPLE} > {output[0]}"
              " && {tools[seqtk]} sample -s {seed} {input[1]} {SUBSAMPLE} > {output[1]}"
              " && [[ -s {output[0]} ]] && [[ -s {output[1]} ]]")
        check_identical_fastq_headers(output[0], output[1])

# Combine the subsampled sequences for each pair
rule combine_each_pair:
    input: expand("temp/01_subsample/{samples}.{{pair}}.{subsample}.fastq",
                  samples=SAMPLES,
                  subsample=SUBSAMPLE)
    output: COMBINED_DATA
    run:
        check_files_arent_empty(input)
        shell("cat {input} > {output} && [[ -s {output} ]]")

# Run fastqc on combined subsampled reads
rule fastqc_each_pair:
    input: COMBINED_DATA
    output: FASTQC_DATA
    threads: MAX_THREADS
    run:
        check_files_arent_empty(input)
        shell("{tools[fastqc]} -t {threads} {input}")
        # Output for fastqc can only be specified to an output directory
        # Therefore, && [[ -s {output} ]] won't work.
        check_files_arent_empty(output)
