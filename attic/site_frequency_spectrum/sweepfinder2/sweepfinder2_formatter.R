#!/usr/bin/env Rscript
############################################################
##
## Author: Joe Colgan           Name: sweepfinder2_format.R
##
## Purpose:
## This script takes an allele frequency file as input,
## rearranges the columns, makes some simple calculations
## and output a tab-delimited text file for use with the
## software program, SweepFinder2.
############################################################
## Take arguments from the command line:
args <- commandArgs(TRUE)
input <- args[1]            # Assign first argument as input
output <- args[2]           # Assign second argument as output

## The input contains eight columns:
## We want to create an output file that contains four columns
## Column 1: Genomic position (SNP position)
## Column 2: Number of derived alleles (alternative)
## Column 3: Total number of samples
## Column 4: Whether the allele is polarised or not.
## Polarised refers to if the allele is known to be derived or
## ancestral. If unpolarised, this column should be '1'

## Read in input file:
data <- read.table(file = input, header = FALSE)

## Correct number of samples:
data$n <- data$V4 / 2
## Calculate the number of derived alleles:
data$x <- round(data$n * data$V8)
## Add status of folded or not:
data$folded <- 1

## Extract columns of interest:
data_to_export <- data[, c(2,10,9,11)]

## Ensure column names are correct:
colnames(data_to_export) <- c("position",
                      "x",
                      "n",
                      "folded")
## Write table to output:
write.table(data_to_export,
            file = output,
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")
