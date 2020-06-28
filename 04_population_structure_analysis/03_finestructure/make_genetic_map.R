#!/usr/bin/env Rscript
##############################################################################
## Author: Joe Colgan                   Program: make_recomb_file.R
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
## Check files ending in *genetic.map.tmp.input.txt

## Provide access to a copy of the command line arguments supplied when R session is invoked
args <- commandArgs(TRUE)
recomb_input <- args[1]           # Assign first argument as input
output <- args[2]
message(paste("will output to", output))

## Update chromosome name:
chromosome <- gsub(pattern = ".recomb.txt",
                   replacement = "",
                   x = recomb_input)

## Assign effective population size
## Effective population size calculated by B. Viera using PSMC
Ne <- 200000

make_genetic_map <- function(){
        ## Read file contents into variable:
        data <- read.table(recomb_input,
                           header = FALSE,
                           col.names = c("locus",
                                         "rho"))
        ## Convert rho values
        data$rho <- data$rho * (100 / ( 3 * Ne ))
        ## Calculate the cumulative values of the rho values
        data$rho <- cumsum(data$rho)
        ## Add snp names
        data$snp <- "."
        data$chrom <- chromosome
        ## Rearrange columns:
        recomb.file <- data[, c(4, 3, 1, 2)]
        ## Write to file
        write.table(recomb.file,
                    file = output,
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE,
                    sep = "\t")
}

## Call the function
make_genetic_map()
