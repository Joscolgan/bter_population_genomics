#!/usr/bin/env R
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: median_counter.R
#
# Date: 29/04/2016
#
##############################################################################

args <- commandArgs(TRUE)
input <- args[1]                                                # Assign first argument as input
output <- args[2]           # Assign second argument as output

message(paste("will output to", output))

## The purpose of this script is to take in a text file containing three columns:
##      - Genomic scaffold - base - read depth
##        NC_89788              1       6
##              ..              ..      ..
## and calculate the median for the third column "read depth"

# Define function 'median_counts'
median_counts<-function() {
        data<-read.table(input, header=FALSE) ## Read in input from command line
        data.colMedians <- median(as.matrix(data[,3]))
        write.table(data.colMedians, file = output, row.names = F, col.names = F, sep="\t")
}

## Call the function
median_counts()
