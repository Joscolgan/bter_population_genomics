#!/usr/bin/env Rscript
##############################################################################
## Author: Joe Colgan                   Program: make_genetic_map.R
##
## Date: 15/03/2015
##
## Purpose:
## To:
## - Read in arguments as {input} and {output}
## - Parse input data for columns containing population-scaled recombination rates (rho)
## - Convert the rho values to centimorgans
## - Generate cumulative sums
## - Add locus identifiers (beneficial for use with selscan EHH)
## - Output to file
##
##############################################################################
## Input files:
## Check files ending in *input_for_genetic_map.txt

## Provide access to a copy of the command line arguments supplied when R session is invoked
args <- commandArgs(TRUE)
input <- args[1] 						# Assign first argument as input
output <- args[2]                                               # Assign second argument as input
message(paste("will output to", output))

## Assign effective population size
## Effective population size calculated by B. Viera using PSMC
Ne <- 200000

## Define function for generating the genetic map
make_genetic_map <- function(){
        ## Read file contents into variable
        data <- read.table(input)
        ## Name columns
        colnames(data) <- c("chrom", "locus", "rho", "physical_distance")
        ## Convert rho values
        ## Conversion taken from here: https://sourceforge.net/p/ldhat/mailman/ldhat-help/thread/072A5A8101E2734E858CEE665023B2986C80287A@SRVUNIMBX01.uni.au.dk/
        data$rho <- data$rho * (100/(3*Ne))
        ## Calculate the cumulative values of the rho values
        data$rho<-cumsum(data$rho)
        ## Add locus names
        data$locus<-1:nrow(data)
	## Add a one to last row physical position to remove duplication:
	data$physical_distance[nrow(data)] <- data$physical_distance[nrow(data)] + 1
        ## Write to file
        write.table(data, file = output,
        row.names = FALSE,
        col.names= FALSE,
        quote = FALSE,
        sep="\t")
}

## Call the function
make_genetic_map()
