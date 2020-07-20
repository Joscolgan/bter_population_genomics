#!/usr/bin/env Rscript
##############################################################################
## Author: Joe Colgan                   Program: make_absolute.R
##
## Date: 27/02/2019
##
## Purpose:
## To:
## - Read in arguments as {input} and {output} from the command line.
## - Input is a BED File containing four fields.
## - Convert a specific field (the fourth) within input containing standardised nsl/ihs values
##   to absolute values.
## - Output to file retaining the same file structure as input.
##
##############################################################################
args <- commandArgs(TRUE)
input <- args[1]                                                # Assign first argument as input
output <- args[2]           # Assign second argument as output

## read in input file and convert a single field to an absolute value:
data <- read.table(file = input, header = FALSE)

## Convert to absolute value:
data$V4 <- abs(data$V4)

## Write to output:
write.table(data,
            file = output,
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t")
