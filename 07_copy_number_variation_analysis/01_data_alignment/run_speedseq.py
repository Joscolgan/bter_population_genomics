#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: run_speedseq.py
#
# Date: 06/05/2016
#
##############################################################################
# Import modules
import os.path
import re

from helper_functions import *

# The purpose of this script is to take two FASTQ files and align against a user-defined
# reference genome to generate an alignment files for:
# - Pairs of mapped reads
# - Split reads
# - Discordant reads
# Final output is three BAM files containing the above information.
#
# To achieve final output, the Snakefile contains custom-defined rules (see below) that
# outline commands to execute sequentially to take custom-defined input(s) and generate
# final output (as defined in rule all).

# For this specific script, rules are defined as follows:
# rule all:
#   - Defines the expected final output of the Snakefile.
# rule align_to_bam:
#   - For each input, align reads against the reference genome using speedseq.

##############################################################################
# Sample information
##############################################################################
# Bumblebee (Bombus terrestris) males were collected summer 2014
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
#   speedseq
#   Adjust the speedseq.config file for paths of software tools required by 
#   speedseq (speedseq will come with packages of software it uses (e.g. bwa, freebaye,
#   lumpy) but the user can download newer/older versions if required and modify
#   the speedseq.config file to direct speedseq to the correct directories containing
#   software of choice. 
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
REFERENCE   = "./GCF_000214255.1_Bter_1.0_genomic_plus_mito_partial_genome.fasta"

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
# Assign path to input reads to align
FILTERED_READ_DATA       = ["../2016-10-21_Bter_read_filtering/temp_test/09_final_clean/{{samples}}.R1.cleaned.fastq",\
                            "../2016-10-21_Bter_read_filtering/temp_test/09_final_clean/{{samples}}.R2.cleaned.fastq"]

# Output bwa-aligned BAM files here
ALIGNED_DATA             = ["./{samples}.bam",\
                            "./{samples}.discordants.bam",\
                            "./{samples}.splitters.bam"]

##############################################################################
# Define binaries in context of path relative to Snakefile
##############################################################################
# binaries
# Align lines of code using 'Assign Align' using cmd+shift+p and selecting 'align'
# Create dictionaries for directories and  tools'
dirs  = {}
dirs['project'] = os.path.abspath('../../../')
dirs['src']     = os.path.join(dirs['project'], 'bin')

# Create an empty dictionary called 'tools'
tools = {}
tools['speedseq']                = os.path.join(dirs['src'], 'speedseq')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################
# First rule is list the final output
rule all:
    input: expand(ALIGNED_DATA, samples=SAMPLES)
            
# Each rule will specify an intermediate step
# Align sample tio the reference genome using bwa-mem
rule align_to_bam:
    input:  forward=expand(FILTERED_READ_DATA[0]),\
            reverse=expand(FILTERED_READ_DATA[1])
    output: ALIGNED_DATA
    run:
        align_list=re.split('/|.R', (''.join(input.forward)))[4]
        print(align_list)
        shell("{tools[speedseq]} align -o {align_list} -R '@RG\tID:{align_list}\tSM:{align_list}\tLB:lib1' -t 10 \
              {REFERENCE} {input.forward} {input.reverse}")

