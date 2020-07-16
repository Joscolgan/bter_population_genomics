#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: run_lumpyexpress.py
#
# Date: 06/05/2016
#
##############################################################################
# Import modules
import os.path

from helper_functions import *

# The purpose of this script is to take two FASTQ files, align against a user-defined
# reference genome and using lumpyexpress, identify structural variants (SVs) and output.
# Final output is a VCF file containing information on SVs.
#
# To achieve final output, the Snakefile contains custom-defined rules (see below) that
# outline commands to execute sequentially to take custom-defined input(s) and generate
# final output (as defined in rule all).

# For this specific script, rules are defined as follows:
# rule all:
#   - Defines the expected final output of the Snakefile.
# rule align_to_bam:
#   - For each input, align reads against the reference genome using BWA-MEM.
# rule add_RG_information:
#   - For each sample, add read group information to input BAM files using Picard.
#   The following information is added:
#   Read Group ID: RGID=sample_ID
#   Read Group Library: RGLB=lib1
#   Read Group PLatform: RGPL=illumina
#   Read Group Platform Unit: RGPU=unit1
#   Read Group sample name: RGSM=sample_ID
# rule filter_discordants:
#   - For each input, filter discordant reads (i.e. PE reads aligning to different co-ordinates)
#   and output in BAM format.
# rule filter_splitters:
#   - For each input, filter split reads (i.e. reads split across two genomic co-ordinates) and
#   output in BAM format.
# rule sort_discordants:
#   - For each input, sort BAM file based on read co-ordinates and output in BAM format.
# rule sort_splitters:
#   - For each input, sort BAM file based on read-co-ordinates and output in BAM format.
# rule run_lumpyexpress:
#   - For each sample, outputs generated from previous rules (aligned BAMs, sorted discordant BAMs
#   and sorted split read BAMs) are used as input to lumpyexpress to calculate the probability
#   for presence of SVs within each sample. Output is in VCF format.

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
#   bwa
#   samtools
#   picard
#   lumpy-sv
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
# Assign path to reference genome fasta file
REFERENCE   = "../GCF_000214255.1_Bter_1.0_genomic.fna"

# Assign the number of threads to be used
THREADS     = 10

##############################################################################
# Assignment of wildcards to be used within rules
##############################################################################
# Open file and read in contents - one sample per line
with open('samples_list.txt') as samples:
    content = samples.readlines()
    SAMPLES = [samples.rstrip('\n') for samples in content]
    print(SAMPLES)

##############################################################################
# Specify all input/output files in terms of sample wildcards
##############################################################################
# Output sorted splitter BAM file here
BAM_DATA     = "{samples}.bam"

# Output sorted discordant BAM file here
DISCORDANT_DATA   = "{samples}.discordants.bam"

# Output unsorted splitter BAM file here
SPLITTER_DATA   = "{samples}.splitters.bam"

# Output the vcf generated by lumpyexpress
LUMPYEXPRESS_OUTPUT      = "combined.all_samples.vcf"

##############################################################################
# Define binaries in context of path relative to Snakefile
##############################################################################
# binaries
# Align lines of code using 'Assign Align' using cmd+shift+p and selecting 'align'
# Create dictionaries for directories and  tools'
dirs  = {}
dirs['project']             = os.path.abspath('../../../../')
dirs['src']                 = os.path.join(dirs['project'], 'src')

# Create an empty dictionary called 'tools'
tools = {}
tools['bwa']                = os.path.join(dirs['src'], 'bwa/bwa')
tools['samtools']           = os.path.join(dirs['src'], 'samtools-1.2/samtools')
#tools['picard']             = os.path.join(dirs['src'], 'picard-tools-1.141/picard.jar')
tools['lumpyexpress']       = os.path.join(dirs['src'], '2016-05-03_lumpySV/lumpy-sv/bin/lumpyexpress')
tools['extractSplitReads']  = os.path.join(dirs['src'], '2016-05-03_lumpySV/lumpy-sv/scripts/extractSplitReads_BwaMem')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################
# First rule is list the final output
rule all:
    input: LUMPYEXPRESS_OUTPUT


# Run lumpyexpress
# Currently an issue with running
rule run_lumpyexpress:
    input: align=expand("{samples}.bam", samples = SAMPLES),\
           split=expand("{samples}.splitters.bam", samples = SAMPLES),\
           discord=expand("{samples}.discordants.bam", samples = SAMPLES)
    output: LUMPYEXPRESS_OUTPUT
    run:
        check_files_arent_empty(input)
        align_list = ",".join(map(str, input.align))
        split_list = ",".join(map(str, input.split))
        discord_list = ",".join(map(str, input.discord))
        shell("{tools[lumpyexpress]} -B {align_list} \
                                     -S {split_list} \
                                     -D {discord_list} \
                                     -o {output} \
                                     -P \
                                     && [[ -s {output} ]]")

