#!/usr/bin/env Rscript
##############################################################################
# Author: Joe Colgan                   Program: ipr_term_conversion.R
#
# Date: 16/08/2017
#
# Purpose:
# This script takes one tab delimited text file as input, which contains two columns:
# first column: Seqname; second column: Interproscan (IPR) domains.
# The script concatenates IPR terms associated to each individual protein sequence (Seqname),
# and generates a nonredundant list of IPR terms for each proteins.
# IPR terms associated with protein were assigned to encoding genes.
#
##############################################################################

## Set working directory:
## Provide access to a copy of the command line arguments supplied when R session is invoked
args      <- commandArgs(TRUE)
input     <- args[1] 	    # Assign first argument as input
output    <- args[2]        # Assign second argument as output

## Load data into dataframe:
ensembl_input <- read.table(input,
                            stringsAsFactors = FALSE)

## Aggregate the list, which concatenates the second column:
mymerge <- function(x) {
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
output_df <- aggregate(ensembl_input[-1],
                       by = list(ensembl_input$V1),
                       mymerge)

## Write to output:
write.table(x = output_df,
            file = output,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE)
