#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: median_CNV_depth.shufflebed.ad.py
#
# Date: 19/04/2016
#
##############################################################################
# Import modules
import os.path

from helper_functions import *

# """
# This script takes an alignment BAM file per sample and intersects with multiple BED files
# containing the genomic positions of individual putative duplication events (generated through
# CNVnator). The intersection process outputs individual BAM files containing reads aligned
# to each putative duplication (bin). The number of aligned reads per BAM are calculated for
# each genomic base per individual. The median number of aligned reads is calculated and output
# to a text file for loadiing into R.

# To achieve final output, the Snakefile contains custom-defined rules (see below) that
# outline commands to execute sequentially to take custom-defined input(s) and generate
# final output (as defined in rule all).

# For this specific script, rules are defined as follows:
# rule all:
#   - Defines the expected final output of the Snakefile.
# rule index_genome:
#    - A single fasta file containing genome of interest is provided as input.
#      Input fasta file is indexed by bowtie2.

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
# Prior to use
##############################################################################
# To run alignment_to_coverage.py:
#  1. Download and install the following software:
#   samtools
#   bedtools2
#
#  2. Ensure helper_functions.py is within the same directory of the Snakefile.
#
#  3. Assign global variables for use within specific rules
#     Please see section below for further information on variable to be assigned
#
#  4. Assignment of wildcards to be used within the rules
#
#  5. Define all the input and outputs of each rule with respect to the above defined wildcards
#
#  6. Each rule willl take assigned input, global variables and generate defined outputs.
#
#  7. Input data should be formatted in the context of defined wildcards: {sample}.{pair}.fastq
#      For example: 2014_Bter_P_D_14_260_head.R1.fastq
#
# 8. Make a text.file called 'sample_list.txt' and put in same directory as Snakfile.
#     Populate 'sample_list.txt' with names of samples to be analysed.
#       For example:
#       2014_Bter_P_D_14_260_head
#       2014_Blap_P_D_21_412_thorax
#       2014_Bpas_A_D_20_391_thorax
#
##############################################################################
# Assign global variables for use in rules (see below)
##############################################################################
# Assign name for list of duplications
DUPLICATION_LIST = "shuffle_bed/shuffle_dup.NW_removedad.txt"

##############################################################################
# Assignment of wildcards to be used within rules
##############################################################################
# Open file and read in contents - one sample per line
with open('samples_list.txt') as samples:
    content = samples.readlines()
    SAMPLES = [samples.rstrip('\n') for samples in content]
    print(SAMPLES)

# Read in the content of the duplication list - one duplication per line
with open('shuffle_bed/shuffle_dup.NW_removedad.txt') as shuffles:
    dup_content = shuffles.readlines()
    SHUFFLES = [shuffles.rstrip('\n') for shuffles in dup_content]
    print(SHUFFLES)

##############################################################################
# Specify all input/output files in terms of sample wildcards
##############################################################################
# Assign path for input alignment BAM files.
ALIGNED_DATA            = "./{samples}.bam"

# Assign path for BED files containing CNV information.
DUPLICATED_DATA         = "shuffle_bed/{shuffles}.bed"

# Output intersected BAM files here.
INTERSECTED_DATA        = "shuffle_median_temp_ad/01_intersect_bam/{shuffles}.{samples}.bam"

# Output depth counts here.
DEPTH_DATA              = "shuffle_median_temp_ad/02_depth_bam/{shuffles}.{samples}.depth.txt"

# Output mean counts per sample here.
CALCULATED_MEDIAN_DATA  = "shuffle_median_temp_ad/03_median_counts/{shuffles}.{samples}.median_depth.txt"

# Output combined mean data here.
COMBINED_MEDIAN_DATA    = "shuffle_median_temp_ad/04_combined_counts/{shuffles}.combined.median_depth.txt"

# Output combined CNV mean data here.
COMBINED_CNV_DATA       = "shuffle_median_temp_ad/05_combined_CNV_means/combined.median_depth.all_CNVs.txt"

# Output final data file here.
FINAL_DATA              = "shuffle_median_temp_ad/06_final/final_combined.median_depth.all_dups.txt"

##############################################################################
# Define binaries in context of path relative to Snakefile
##############################################################################
# binaries
# Align lines of code using 'Assign Align' using cmd+shift+p and selecting 'align'
# Create dictionaries for directories and  tools'
dirs  = {}
dirs['project']        = os.path.abspath('../../../../')
dirs['src']            = os.path.join(dirs['project'], 'src')

# Create an empty dictionary called 'tools'
tools = {}
tools['intersect']     = os.path.join(dirs['src'], 'bedtools2/bin/intersectBed')
tools['samtools']      = os.path.join(dirs['src'], 'samtools-1.2/samtools')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################
# First rule is list the final output
rule all:
    input: expand(FINAL_DATA)

# Each rule will specify an intermediate step
# Intersect BAM file using BED files populated with CNV genomic co-ordinates
rule intersect_bam:
    input:  ALIGNED_DATA, DUPLICATED_DATA
    output: INTERSECTED_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[intersect]} -abam {input[0]} -b {input[1]} > {output} && [[ -s {output} ]]")

# Count depth of coverage for reads aligned within intersected BAM files.
rule count_depth:
    input:  INTERSECTED_DATA
    output: DEPTH_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[samtools]} depth {input} > {output} && [[ -s {output} ]]")

# Calculate the median of aligned reads for each putative CNV.
rule calculate_median:
    input:  DEPTH_DATA
    output: CALCULATED_MEDIAN_DATA
    run:
        check_files_arent_empty(input)
        shell("Rscript median_counter_test.R {input} {output} && [[ -s {output} ]]")

## Original code used to calculate the mean depth.
## shell("awk '{{ sum += $3 }} END {{ print sum/NR}}' {input} > {output} && [[ -s {output} ]]")

# Paste the means for each sample per row
rule paste_means:
    input:  expand("shuffle_median_temp_ad/03_median_counts/{{shuffles}}.{samples}.median_depth.txt", samples = SAMPLES)
    output: COMBINED_MEDIAN_DATA
    run:
        check_files_arent_empty(input)
        shell("paste {input} > {output} && [[ -s {output} ]]")

# Combine the CNV mean rows
rule combine_CNVs:
    input:  expand("shuffle_median_temp_ad/04_combined_counts/{shuffles}.combined.median_depth.txt", shuffles = SHUFFLES)
    output: COMBINED_CNV_DATA
    run:
        check_files_arent_empty(input)
        shell("cat {input} > {output} && [[ -s {output} ]]")

# Combine counts with mean row name
rule add_CNV_names:
    input:  COMBINED_CNV_DATA
    output: FINAL_DATA
    run:
        check_files_arent_empty(input)
        shell("paste {DUPLICATION_LIST} {input} > {output} && [[ -s {output} ]]")
