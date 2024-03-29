#!/usr/bin/env Rscript
#######################################################################
##
## Author: Joe Colgan           Program: plot_pca.R
##
## Date:   2019/05/10
##
## Purpose:  
## This script takes a VCF file as input (VCF generated by Freebayes) and using gdsfmt generates a 
## genotype matrix.
## Using this genotypic matrix, SNPRelate performs LD pruning, PCA and isolation by descent to
## assess relatedness amongst individuals.
## The output is a number of ggplots for:
## - A scatterplot of first two principal components.
## - A scatterplot for IBD co-efficients for each individual in the dataset.
## These plots are output as .PNG files.
#######################################################################

# Load libraries; install from scratch if needed
libraries <- c("gdsfmt",
               "SNPRelate",
               "ggplot2",
               "ggpubr",
               "lintr",
               "reshape2")
for (lib in libraries) {
    if (require(package = lib, character.only = TRUE)) {
        print("Successful")
    } else {
        print("Installing")
        source("https://bioconductor.org/biocLite.R")
        library(lib, character.only = TRUE )
    }
}

## Create output directory:
dir.create("results")

## Take arguments from the command line:
## Provide access to a copy of the command line arguments supplied when R session is invoked
args <- commandArgs(TRUE)
input <- args[1] 						# Assign first argument as input
output <- args[2]						# Assign second argument as input
message(paste("will output to", output))

## Load data and perform PCA:
## Plot PCA of larger number of samples (n=46):
## Assign the input name:  
input_vcf <- input

## Convert VCF to genome data structure (gds) file format
## gds is a container for storing annotation data and SNP genotypes.
## In this format, each byte encodes up to four SNP genotypes thereby reducing file size and access time.
## The GDS format supports data blocking so that only the subset of data that is being processed needs
## to reside in memory. GDS formatted data is also designed for efficient data access to large datasets.

## Use the name of input to designate the name for the GDS file format.
input_vcf_gds <- gsub(".vcf", ".gds", input_vcf)

## Reformat VCF into GDS format
snpgdsVCF2GDS(input_vcf,
              input_vcf_gds,
              method = "biallelic.only")

## Extract genotype information:
input_vcf_gds_genofile <- snpgdsOpen(input_vcf_gds)

## Extract genotype information and generate genotypic matrix:  
input_vcf_gds_genofile_matrix <- read.gdsn(index.gdsn(input_vcf_gds_genofile,
                                                      "genotype"))
                                                      
## For performing PCA, it is recommended to perform LD pruning:
snpset <- snpgdsLDpruning(input_vcf_gds_genofile,
                          ld.threshold = 0.1,
                          autosome.only = FALSE)

## Get all selected SNP ids:
snpset.id <- unlist(unname(snpset))

## Perform PCA and store output in variable:  
pca        <- snpgdsPCA(input_vcf_gds_genofile,
                          num.thread = 20,
                          snp.id = snpset.id,
                          autosome.only = FALSE)

## Calculate percentages for principal components:  
pc_percent      <- pca$varprop * 100
## Print rounded up percentages:
print(round(pc_percent, 2))

## Extract sample names from input:  
sample_id     <- read.gdsn(index.gdsn(input_vcf_gds_genofile,
                                      "sample.id"))

## Read in file containing treatment information. 
## This file should contain 'treatment information' for each sample and be in the same order as each sample. 
## i.e. the first treatment information would be assigned to the first sample, etc., etc.
## Read in table:
pop_code_df <- read.table(file = "input/population_information.txt",
                          header = FALSE)

## Rename columns:
colnames(pop_code_df) <- c("sample",
                           "site_number",
                           "landuse",
                           "original_site_number",
                           "latitude")

## 
pop_code <- pop_code_df$landuse

## Generate a dataframe containing a column for:
## 1) Sample name
## 2) Treatment condition
## 3) First principal component of interest
## 4) Second principal component of interest
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample_id)],
                  EV1 = pca$eigenvect[, 1],    # the first eigenvector
                  EV2 = pca$eigenvect[, 2],    # the second eigenvector
                  stringsAsFactors = FALSE)

## Rename the samples for plotting:  
## Renaming the sample.ids
new_names <- list()

## Create a list:
for (name in sample_id){
        new_names[name] <- paste(strsplit(name, "_")[[1]][5],
        "_",strsplit(name, "_")[[1]][6],
        sep = "")
}

## Unlist as a character string and update sample ids:
tab$sample.id <- as.character(unlist(new_names))

## Add another column for site_id:
tab$site_id <- pop_code_df$site_number

## Add column for latitude:
tab$latitude <- pop_code_df$latitude

## Check status of whether there are two or one individual per site:
new_site_id <- vector()

for (i in 1:length(tab$site_id)){
        status <- duplicated(tab$site_id)[i]
        if (status==FALSE){
                print("Not duplicated")
                new_site_id <- c(new_site_id, paste(tab$site_id[i], "A", sep=""))
        }
        else {
                print("Duplicated")
                new_site_id <- c(new_site_id, paste(tab$site_id[i], "B", sep=""))
        }
}

## Update site information: 
tab$site_id <- new_site_id

## Alternatively, to plot with just sample ids and not points:
names_plot<- ggplot(tab,
                    aes(x = tab$EV1,
                        y = tab$EV2,
                        #color = tab$pop,
                        fill = tab$latitude)) +
        #geom_point(shape = 16, size = 4, alpha = 0.4, show.legend = T) +
        theme_minimal() +
        xlab(paste("Principal component", 1," ","(",round(pc_percent[1], 2), "%", ")", sep = "")) +
        ylab(paste("Principal component", 2," ","(",round(pc_percent[2], 2), "%", ")", sep = "")) + 
        geom_text(size=6, position=position_jitter(width = 0.04, height = 0.04),
        aes(label = tab$site_id, color = tab$latitude)) +
        scale_fill_distiller(palette = "Blues")
        #scale_color_manual(values = c("blue", "orange", "red"))

## Adjust the size of the labels:
names_plot <- names_plot +
              theme(axis.text=element_text(size = 20),
                    axis.title=element_text(size = 20,
                                            face = "bold")) +
                    theme(legend.position = "none")

## Plot:
ggsave(file = output,
       dpi = 600,
       width = 10,
       height = 10)
