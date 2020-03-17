#!/usr/bin/env Rscript
##############################################################################
# Author: Joe Colgan                   Program: GO_term_conversion.R
#
# Date: 16/08/2017
#
# Purpose:
# This script takes one tab delimited text file as input, which contains two columns:
# first column: Seqname; second column: GOs.
# The script concatenates GO terms associated to each individual protein sequence (Seqname),
# and generates a nonredundant list of GO terms for each proteins.
# GO terms associated with protein were assigned to encoding genes.
#
##############################################################################

## Take arguments from the command line:
args <- commandArgs(TRUE)
input <- args[1] 						# Assign first argument as input
output <- args[2]           # Assign second argument as output

## Load data into dataframe:
input <- read.table(input,
stringsAsFactors = FALSE)

## Aggregate the list, which concatenates the second column:
mymerge <- function(x){
        all_in_one <- paste(unlist(x),
        sep = ",",
        collapse = ",")
        split_term <- unlist(strsplit(all_in_one,
        split = ","))
        return(paste(unique(split_term),
        sep = ",",
        collapse = ","))
}

## Run function:
output_file <- aggregate(input[-1],
              by = list(input$V1), mymerge)

## Write output to file:
write.table(output_file,
            file = output,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
