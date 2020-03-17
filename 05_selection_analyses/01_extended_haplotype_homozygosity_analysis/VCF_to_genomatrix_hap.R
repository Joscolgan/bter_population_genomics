#!/usr/bin/env Rscript
##############################################################################
## Author: Joe Colgan                   Name: VCF_to_genomatrix_hap.R
## 
## Date: 2017-01-24
## 
## Purpose:
## Code to take VCF as input and generate genotypic matrix for each.
## Then to output genotype matrices
##
##=================================================================================

# Load libraries; install from scratch if needed
libraries <- c("gdsfmt", "SNPRelate", "MASS",
               "corrplot", "ggplot2", "reshape2")
for (lib in libraries) {
    if (require(package = lib, character.only = TRUE)) {
        print("Successful")
    } else {
        print("Installing")
        source("https://bioconductor.org/biocLite.R")
        biocLite(pkgs = lib)
        library(lib, character.only = TRUE )
    }
}

# Load libraries
#library(gdsfmt, lib='/data/home/btw928/R/x86_64-unknown-linux-gnu-library/3.2/')
#library(SNPRelate, lib='/data/home/btw928/R/x86_64-unknown-linux-gnu-library/3.2/')
#library(MASS, lib='/data/home/btw928/R/x86_64-pc-linux-gnu-library/3.3/')
#library(corrplot, lib='/data/home/btw928/R/x86_64-unknown-linux-gnu-library/3.2/')
#library(ggplot2, lib='/data/home/btw928/R/x86_64-unknown-linux-gnu-library/3.2/')
#library(reshape, lib='/data/home/btw928/R/x86_64-unknown-linux-gnu-library/3.2/')

## Provide access to a copy of the command line arguments supplied when R session is invoked
args <- commandArgs(TRUE)
input <- args[1]            # Assign first argument as input
input.gds <- args[2]            # Assign second argument as input.gds file
output_hap <- args[3]           # Assign third argument as hap output
output_map <- args[4]           # Assign fourth argument as map output

message(paste("will output to", output_hap))


# Define function 'vcf_to_genomatrix'
vcf_to_genomatrix_hap<-function(data,
                                data.gds,
                                output_hap,
                                output_map) {
        # Reformat the input
        snpgdsVCF2GDS(data, data.gds, method= "biallelic.only")
        ## Generate a summary of SNP information
        print(snpgdsSummary(data.gds))
        ## Open the GDS file
        genofile <- snpgdsOpen(data.gds)
        ## Output "genotype" matrix and store within variable 'g'
        genotype_matrix <- read.gdsn(index.gdsn(genofile, "genotype"))
        ## Format change '2' to '1'
        genotype_matrix[genotype_matrix==2] <- 1
        snp_chromosome <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
        snp_position <- read.gdsn(index.gdsn(genofile, "snp.position"))
        output_df <- as.data.frame(cbind(snp_chromosome, snp_position))
        output_df$snp_name <- "."
        output_df$snp_position_2 <- snp_position
        output_df <- output_df[ ,c(1,3,2,4)]
        ## Write matrix to file
        write.matrix(genotype_matrix, file = output_hap, sep = " ")
        ## SNP positions to file
        write.table(output_df, file = output_map, sep = "\t", col.names = FALSE,
                    row.names = FALSE,
                    quote = FALSE)
}

vcf_to_genomatrix_hap(data = input,
                      data.gds = input.gds,
                      output_hap = output_hap,
                      output_map = output_map)
