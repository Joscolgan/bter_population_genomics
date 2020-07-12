#!/usr/bin/env Rscript
##############################################################################
# Author: Joe Colgan                   Program: GO_term_conversion.R
#
# Date: 16/08/2017
#
# Purpose:
# This script takes one tab delimited text file as input, which contains two columns:
# first column: Seqname; second column: GOs.
# The script concatenates GO terms associated to each individual locus name,
# and generates a nonredundant list of GO terms for each locus.
#
##############################################################################

## Take input from the command line:
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

input <- args[1]
output <- args[2]

## Load data into dataframe:
go_term_input <- read.table(file = input,
stringsAsFactors = FALSE)

## Aggregate the list, which concatenates the second column:
mymerge <- function(x) {
        all_in_one <- paste(unlist(x), sep=",", collapse=",")
        split_term <- unlist(strsplit(all_in_one, split=","))
        return(paste(unique(split_term), sep=",", collapse=","))
}

## Run function:
aggregated_terms <- aggregate(go_term_input[-1], by=list(go_term_input$V1), mymerge)

write.table(aggregated_terms,
	file = output,
	sep = "\t",
	quote = FALSE,
	row.names = FALSE,
	col.names = FALSE)
