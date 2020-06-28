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
positions <- args[1]
recomb_input <- args[2]           # Assign first argument as input
output <- args[3]
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
        #data$rho <- cumsum(data$rho)
        ## Duplicate last row:
        data <- rbind(data, data[nrow(data),])
        ## Read in positions and update:
        map <- read.table(positions,
                          header = FALSE,
                          col.names = c("chrom",
                                        "snp",
                                        "recomb_zero",
                                        "locus"))
        ## Add snp names
        #data$snp <- 1:nrow(data)
        #data$chrom <- chromosome
        data$locus <- map$locus
        ## Convert final recombination value to 0:
        data$rho[nrow(data)] <- 0
        colnames(data) <- c("start.pos",
                            "recom.rate.perbp")
        ## Write to file
        write.table(data,
                    file = output,
                    row.names = FALSE,
                    col.names = TRUE,
                    quote = FALSE,
                    sep = "\t")
}

## Call the function
make_genetic_map()
