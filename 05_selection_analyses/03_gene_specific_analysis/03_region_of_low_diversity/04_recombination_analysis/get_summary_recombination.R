#!/usr/bin/env Rscript
#######################################################################
##
## Author: Joe Colgan           Program: get_summary_recombination.R
##
## Date:   2018/07/03
##
## Purpose:
## This script calculates the summary statistics for population scaled 
## recombination rate calculated for each SNP that overlaps with a 100 kb 
## genomic region of significant low diversity.
##
#######################################################################

## Provide access to a copy of the command line arguments supplied when R session is invoked
args <- commandArgs(TRUE)
input <- args[1] 						# Assign first argument as input
output <- args[2]						# Assign second argument as input
message(paste("will output to", output))

## Calculate summary statistics for 
## Update column names:
input_df <- read.table(file=input, header=FALSE)

colnames(input_df) <- c("chrom", "start", "end", "rho")

## Calculate summary statistics for rho (population-scaled recombination) values:
summary_stats <- data.frame()
summary_stats <- as.data.frame(rbind(as.character(unlist(summary(input_df$rho)))))

## Update column names:
colnames(summary_stats) <- c("min.",
                             "1st_Q",
                             "Median",
                             "Mean",
                             "3rd_Q",
                             "Max")

## Reduction from genome-wide rho mean ("0.000349814")
summary_stats$genome_wide_mean <- 0.000349814

summary_stats$diff_from_mean <-  summary_stats$genome_wide_mean/as.numeric(as.character(unlist(summary_stats$Mean)))

## Write summary statistics to output:
write.table(x=summary_stats,
            file=output,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t",
            quote = FALSE)
