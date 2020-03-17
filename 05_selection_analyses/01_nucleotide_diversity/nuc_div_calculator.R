#!/usr/bin/env Rscript
##########################################################################################
##
## Author: Joe Colgan                           Program: nuc_div_calculator.R
##
## Date: 04-03-2018
##
## Purpose:
##  - Read in text file containing information of VCF file for each chromosome, as well
##    as start base position and end base position for each corresponding chromosome.
##  - Calculate diversity and neutrality statistics for each chromosome.
##  - Output a table containinig -
##      a) Chromosome
##      b) start position of genomic window
##      c) end position of genomic window
##      d) nucleotide diversity for user-defined genomic window
##      e) number of segregating sites for user-defined genomic window
##      f) Tajima's D for use-defined genomic window
##      g) midpoint position of genomic window
##
##########################################################################################

input <- "./chromosome_info.txt"
bin_width <- 100000
input_gff <- "../gff/Bombus_terrestris.Bter_1.0.44.gff"
output <- "./bter_all_chrom_100kb_approx_F_jump_size_zero.txt"

# Load libraries; install from scratch if needed
libraries <- c("ggplot2",
"PopGenome",
"stringr",
"lintr")
for (lib in libraries) {
    if (require(package = lib, character.only = TRUE)) {
        print("Successful")
    } else {
        print("Installing")
        source("https://bioconductor.org/biocLite.R")
        library(lib, character.only = TRUE )
    }
}

#install.packages("WhopGenome")

## Load data
## Read in the tabix indexed VCF
## The input parameters include:
## 1) Tabix-indexed bgzipped file
## 2) numcols: number of SNPs that should be read in as a chunk
## 3) tid: The ID of the chromosome e.g. "NC_015764.1"
## 4) Frompos: From position e.g. 1
## 5) topos: End positon (probably length of chromosome)
## An additional parameter is to include a gff

## Scan in information about chromosomes:
chrom_data <- read.table(input, header=F)
chrom_data <- as.matrix(chrom_data)

## Divide each chromosome into sliding windows (window size: 10kb, jump size:5kb)
##
count <- 0
nucdiv_combined_df <- data.frame()

## Define jumpwidth:
jump_width <- bin_width
print(jump_width)
print(chrom_data)

## For loop, cycle through each chromosomes and calculate nucleotide diverity within
## genomic windows of user-defined size:
for (item in 1:nrow(chrom_data)){
        ## Extract name of chromosome:
        chrom <- gsub(".recode.vcf.bgz", "", chrom_data[item,][1])
        ## Read in VCF file:
        chrom_data_item <- readVCF(chrom_data[item,][1], bin_width, chrom, chrom_data[item,][2], chrom_data[item,][3], gffpath = input_gff, approx = FALSE, include.unknown=TRUE)
        ##
        chrom_data_item      <- sliding.window.transform(chrom_data_item, bin_width, jump_width, type=2, whole.data = F)
        ## Extract region names, which contains window coordinates:
        region_names         <- chrom_data_item@region.names
        region_names.df      <- str_split_fixed(region_names, " ", 4)[,4]
        region_names.df      <- gsub(" :", "", region_names.df)
        region_names.df      <- str_split_fixed(region_names.df, " - ", 2)
        ## Calculation of neutrality statistics:
        chrom_data_item      <- neutrality.stats(chrom_data_item)
        ## Calculation of diversity statistics:
        chrom_data_item      <- diversity.stats(chrom_data_item, pi=T)
        ## Calculation of nucleotide diversity:
        nucdiv_item          <- chrom_data_item@nuc.diversity.within
        nucdiv_item          <- nucdiv_item/bin_width
        nucdiv_item.df <- as.data.frame(cbind(nucdiv_item, region_names.df))
        nucdiv_item.df$chrom <- chrom
        nucdiv_item.df$tajima_d <- chrom_data_item@Tajima.D
        ## Subset information on the number of segregating sites:
        nucdiv_item.df$seg_sites <- chrom_data_item@n.segregating.sites
        ## Combine nucleotide diversity statistics across windows:
        nucdiv_combined_df <- rbind(nucdiv_combined_df, nucdiv_item.df)
}

## Rename columns:
colnames(nucdiv_combined_df) <- c("nuc_diversity",
                                  "start",
                                  "end",
                                  "chrom",
                                  "tajima_d",
                                  "seg_sites")

## Rearrange order:
nucdiv_combined_df <- nucdiv_combined_df[c(4,2,3,1,5,6)]

## Calcule median nucletide diversity
nucdiv_median <- median(as.numeric(as.character(nucdiv_combined_df$nuc_diversity)))

## Convert to numeric values:
nucdiv_combined_df$nuc_diversity <- as.numeric(as.character(nucdiv_combined_df$nuc_diversity))
nucdiv_combined_df$start         <- as.numeric(as.character(nucdiv_combined_df$start))
nucdiv_combined_df$end           <- as.numeric(as.character(nucdiv_combined_df$end))
nucdiv_combined_df$seg_sites     <- as.numeric(as.character(nucdiv_combined_df$seg_sites))
nucdiv_combined_df$tajima_d      <- as.numeric(as.character(nucdiv_combined_df$tajima_d))
nucdiv_combined_df$midpoint      <- round((nucdiv_combined_df$start + nucdiv_combined_df$end)/2)

## Export table:
write.table(nucdiv_combined_df,
output,
row.names = FALSE,
col.names = TRUE,
quote = FALSE,
sep="\t")
