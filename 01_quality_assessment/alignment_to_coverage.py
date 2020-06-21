#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: alignment_to_coverage.py
#
# Date: 02/06/2015
#
##############################################################################
# Import modules
import os.path

from helper_functions import *

# This script takes two fastq files (i.e. pairs) per sample(s) as input and aligns
#  input reads against indexed reference genome. It then calculates, parses and plots
#  percantage of genomic regions with low coverage (5x > per sequence base).
# Final output is a stacked bar plot plotting:
#    (X-axis): Sample name
#    (Y-axis): Percentage coverage of genomic regions

# To achieve final output, the Snakefile contains custom-defined rules (see below) that
# outline commands to execute sequentially to take custom-defined input(s) and generate
# final output (as defined in rule all).

# For this specific script, rules are defined as follows:
# rule all:
#   - Defines the expected final output of the Snakefile.
# rule index_genome:
#    - A single fasta file containing genome of interest is provided as input.
#      Input fasta file is indexed by bowtie2.
# rule align_to_genome:
#   - For each sample, sequences are aligned in pairs to the bowtie2-indexed genome.
#    The rule checks if sequence headers for each input pair contain the same header information.
#    If they differ, an error is raised.
# rule sort_sam_to_bam:
#    - For each sample, aligned SAM file is converted to BAM file and sorted.
# rule calculate_coverage:
#    - For each sample, calculates read depth per base for each genomic scaffold using an input file
#      containing genomic co-ordinates for each genomic scaffold.
#    - In addition, summarises percentage of read depth across entire genome.
# rule parse_low_coverage_regions:
#    - For each input file, parses rows where the second column contains a value less than 5.
#    - Subsequently, parses rows containing 'genome' (which is summary information) and outputs.
# rule reformat_plot_data:
#    - Using a custom R script, reformat the parsed data into a format that can be plotted.
# rule combine_plot_input:
#    - Combine reformatted data for each sample.
# rule plot_stack_charts:
#    - Using the combined reformatted data for each sample, plot a stacked bar plot.

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
#   bowtie2
#   samtools
#   bedtools2
#
#  2. Ensure helper_functions.py is within the same directory of the Snakefile.
# For this particular snakefile, ensure the location of the bash script(genomecov_output_parser.sh)
#  and Rscripts (reformat_for_plot_test.R and plot_stacks.R) are defined.
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

##############################################################################
# Assign global variables for use in rules (see below)
##############################################################################
# Specify the reference genome file to be indexed
GENOME        = "GCF_000214255.1_Bter_1.0_genomic.fna"

# Specify the custom index of Bowtie2-build to use
INDEX         = "Bter_genome"

# Specify the number of threads
MAX_THREADS   = 10

# Specify coverage file for read depth counts
# genomeCoverageBed requires a genome coverage file which is a tab-delimited file
# with two columns: column one: chromosome name; column two: chromosome length
# The file can be generated by running 'samtools faidx <genome.fasta>', which will
# generate an indexed fasta file with the first column containing the chromosome
# name and the second column containing chromosome lenth. The two columns can be
# cut (cut -f 1,2 <file>) into a genome index file like the one below.
COVERAGE_FILE = "./GCF_000214255.1_Bter_1.0_genomic.index.txt"

##############################################################################
# Assignment of wildcards to be used within rules
##############################################################################
# Open file and read in contents - one sample per line
with open('samples_list.txt') as samples:
    content = samples.readlines()
    SAMPLES = [samples.rstrip('\n') for samples in content]
    print(SAMPLES)

# Assign pair information
PAIR          = ["R1", "R2"]

##############################################################################
# Specify all input/output files in terms of sample wildcards
##############################################################################
# Assign output data for indexed genome data
INDEXED_DATA      = ["{index}.1.bt2",\
                     "{index}.2.bt2",\
                     "{index}.3.bt2",\
                     "{index}.4.bt2",\
                     "{index}.rev.1.bt2",\
                     "{index}.rev.2.bt2"]

# Assign path for clean read data for aligning
CLEANED_DATA      = ["{samples}.R1.cleaned.fastq",
                     "{samples}.R2.cleaned.fastq",
                     "{samples}.unpaired.fastq"]

# Output aligned data and unmapped data in these directories, respectively
ALIGNED_DATA      = "filtered_temp/01_aligned/{samples}.sam"

# Output unmapped data here
UNMAPPED_DATA     = "filtered_temp/01_aligned/{samples}_unmapped/"

# Output sorted and merged BAM files here
SORTED_DATA       = "filtered_temp/02_sorted/{samples}.bam"

RG_ANNOTATED_DATA = "filtered_temp/02_sorted/{samples}_RGadded.bam"

# Output text file containing read depth information here
COVERAGE_DATA     = "filtered_temp/03_count/{samples}.read_depth.txt"

# Output text file containing read depth information here
PARSEDPLOT_DATA   = "filtered_temp/03_count/{samples}.read_depth.filtered.txt"

# Output reformatted data here
REFORMED_DATA     = "filtered_temp/04_reform/{samples}.read_depth.reformed.txt"

# Output combined reformatted plot data here
COMBINEDPLOT_DATA = "filtered_temp/04_reform/plot_data.read_depth.combined.txt"

# Output stacked bar plot of coverage assessment here
STACKPLOT_DATA   = "2016-04-17_plotted_data.read_depth.combined-filteredk5.pdf"

##############################################################################
# Define binaries in context of path relative to Snakefile
##############################################################################
# binaries
# Align lines of code using 'Assign Align' using cmd+shift+p and selecting 'align'
# Create dictionaries for directories and  tools'
dirs  = {}
dirs['project'] = os.path.abspath('../../../')
dirs['src']     = os.path.join(dirs['project'], 'src')

# Create an empty dictionary called 'tools'
tools = {}
tools['build']       = os.path.join(dirs['src'], 'bowtie2-2.2.5/bowtie2-build')
tools['bowtie2']     = os.path.join(dirs['src'], 'bowtie2-2.2.5/bowtie2')
tools['samtools']    = os.path.join(dirs['src'], 'samtools-1.2/samtools')
tools['coverage']    = os.path.join(dirs['src'], 'bedtools2/bin/genomeCoverageBed')
tools['picard']      = os.path.join(dirs['src'], 'picard-tools-1.141/picard.jar')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################
# First rule is list the final output
rule all:
    input: expand(STACKPLOT_DATA)

# Each rule will specify an intermediate step
# Index input genome into format applicable for alignment
rule index_genome:
    input:  GENOME
    output: INDEXED_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[build]} {input} {INDEX}")

# Align cleaned sequences against indexed genome
rule align_to_genome:
    input:  expand(INDEXED_DATA, index=INDEX), CLEANED_DATA,
    output: ALIGNED_DATA, UNMAPPED_DATA
    threads: MAX_THREADS
    run:
        align_list=re.split('/|.R', (''.join(output[0])))[2]
        check_files_arent_empty(input)
        shell("{tools[bowtie2]} --rg-id {align_list} \
                                --rg SM:{align_list} \
                                --rg LB:library1 \
                                --rg PL:ILLUMINA \
                                --rg DS:HiSeq2500 \
                                --local -p {threads} \
                                -X 1000 -x {INDEX} \
                                -1 {input[6]} \
                                -2 {input[7]} \
                                -U {input[8]} \
                                -S {output[0]} && [[ -s {output[0]} ]]")


# Convert SAM file to BAM format and sort
rule sort_sam_to_bam:
    input:  ALIGNED_DATA
    output: SORTED_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[samtools]} view -bS {input} \
             | {tools[samtools]} sort - -f {output} && [[ -s {output} ]]")

# Calculate read depth per genomic scaffold base and coverage percentage across individual genomic
# scaffolds
rule calculate_coverage:
    input:  SORTED_DATA
    output: COVERAGE_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[coverage]} -ibam {input} -g {COVERAGE_FILE} > {output} && [[ -s {output} ]]")

# Parse read depth information for each sample and output data for plotting
rule parse_low_coverage_regions:
    input:  COVERAGE_DATA
    output: PARSEDPLOT_DATA
    run:
        check_files_arent_empty(input)
        shell("bash genomecov_output_parser.sh {input} {output} && [[ -s {output} ]]")

# Reformat data in format applicable with stacked bar plots
rule reformat_plot_data:
    input:  PARSEDPLOT_DATA
    output: REFORMED_DATA
    run:
        shell("Rscript reformat_for_plot_test.R {input} {output} && [[ -s {output} ]]")

# Combine plot data for a comparative plots for samples
rule combine_plot_input:
    input:  expand(REFORMED_DATA, samples=SAMPLES)
    output: COMBINEDPLOT_DATA
    run:
        check_files_arent_empty(input)
        shell("cat {input} > {output} && [[ -s {output} ]]")

# Plot combined read depth information
rule plot_stack_charts:
    input:  COMBINEDPLOT_DATA
    output: STACKPLOT_DATA
    run:
        check_files_arent_empty(input)
        shell("Rscript plot_stacks.R {input} {output} && [[ -s {output} ]]")
