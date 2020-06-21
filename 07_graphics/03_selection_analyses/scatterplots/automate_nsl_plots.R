#!/usr/bin/env Rscript
## Load libraries:
libraries <- c("ggplot2",
               "SNPRelate",
               "reshape2",
               "ggpubr",
               "stringr",
               "gggenes")
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

## Take arguments from the command line:
## Provide access to a copy of the command line arguments supplied when R session is invoked
args <- commandArgs(TRUE)
input <- args[1] 						# Assign first argument as input
output <- args[2]

chrom <- gsub("-.*", "", input)
## Print to console:
print(chrom)
start <- as.numeric(strsplit(input, "-")[[1]][2])
print(start)
end <- as.numeric(strsplit(input, "-")[[1]][3])
print(end)

## Define thresholds for nsl and ihs:
nsl_threshold <- 2.56
ihs_threshold <- 2.89

## This script plots nsl score, ihs score, population-scaled recombination rate,
## nucleotide diversity and Tajima's D estimates.

## Read in data:
nsl_input <- paste("./input/nsl_output/",
                   chrom,
                   ".nsl.out.100bins.norm.abs.bed.tmp",
                   sep = "")

ihs_input <- paste("./input/ihs_output/",
                   chrom,
                   ".ihs.out.100bins.norm.abs.bed.tmp",
                   sep = "")

recomb_input <- paste("./input/recombination_rate/,
                     chrom,
                     ".recode.vcf.ldhelmet.ldhelmet.post.txt.tmp.genetic_map.txt.recomb.txt",
                     sep = "")
nsl_data <- read.table(file = nsl_input,
                       col.names = c("chrom",
                                     "position",
                                     "position",
                                     "nsl_score"))
ihs_data <- read.table(file = ihs_input,
                       col.names = c("chrom",
                                     "position",
                                     "position",
                                     "ihs_score"))
recomb_data <- read.table(file = recomb_input,
                          header = TRUE)
## Both nucleotide diversity and Tajima D are in the same input file:
nuc_div_data <- read.table(file = "./input/nucleotide_diversity_output/bter_all_chrom_10kb.txt",
                           header = TRUE)
## Assign the coordinates of the gene to variables plus 100kb
## up and downstream of the gene:
input_data <- read.table(file = input,
                         header = FALSE)

## Subset nsl values for plotting:
nsl_input_subset <- subset(nsl_data,
                            position > start &
                            position < end)
## subset sites above 2.56 to highlight on plot:
nsl_snps_highlight <- subset(nsl_input_subset,
                            nsl_score > nsl_threshold)
## Subset ihs values for plotting:
ihs_input_subset <- subset(ihs_data,
                            position > start &
                            position < end)
ihs_snps_highlight <- subset(ihs_input_subset,
                            ihs_score > 2.89)

## Generate plot of EHH:
nsl_plot <- ggplot(nsl_input_subset,
                   aes(x = position,
                       y = nsl_score)) +
                xlab("Genomic coordinates (bp)") +
                ylab("|nsl| value") +
                geom_rect(aes(xmin = 764617,
                              xmax = 879060,
                              ymin = -Inf,
                              ymax = Inf),
                          fill = "light blue",
                          alpha = 0.008) +
        geom_hline(yintercept = c(2.56),
                   linetype = "dashed",
                   colour = "orange") +
        scale_x_continuous(labels = scales::comma,
                           limits = c(start, end)) +
                geom_point(colour = "black",
                           size = 3,
                           alpha = 0.5) +
                geom_point(data = nsl_snps_highlight,
                           aes(x = position,
                               y = nsl_score),
                           colour = "blue",
                           size = 3,
                           alpha = 0.5) +
                theme_bw() +
        theme(axis.text = element_text(size = 15,
                                       face = "bold"),
              axis.title.y = element_text(size = 15,
                                          face = "bold"),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

## Generate the ihs plot:

ihs_plot <- ggplot(ihs_input_subset,
                   aes(x = position,
                       y = ihs_score)) +
                xlab("Genomic coordinates (bp)") +
                ylab("|ihs| value") +
                geom_rect(aes(xmin = 764617,
                              xmax = 879060,
                              ymin = -Inf,
                              ymax = Inf),
                          fill = "light blue",
                          alpha = 0.008) +
        geom_hline(yintercept = c(2.89),
                   linetype = "dashed",
                   colour = "orange") +
        scale_x_continuous(labels = scales::comma,
                           limits = c(start, end)) +
                geom_point(colour = "black",
                           size = 3,
                           alpha = 0.5) +
                geom_point(data = ihs_snps_highlight,
                           aes(x = position,
                               y = ihs_score),
                           colour = "blue",
                           size = 3,
                           alpha = 0.5) +
                theme_bw() +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 10,
                                        face = "bold"))

## Plot recombination rate for the region:

## Subset region of genome to plot:
recomb_data_subset <- subset(recomb_data,
                            start.pos > start &
                            start.pos < end)
recomb_plot <- ggplot(recomb_data_subset,
                      aes(x = start.pos,
                          y = recomb.rate.perbp)) +
               geom_point()

## Subset genomic region to plot nucleotide diversity:
## Calculate median nucleotide diversity:
nucdiv_median <- median(nuc_div_data$nuc_diversity)

## Subset corrected genomic regions to plot:
nuc_div_subset <- subset(nuc_div_data,
                         chrom == chrom &
                         start > start &
                         end < end)
options(scipen = 999)
## Plot:
nuc_div_plot <- ggplot(nuc_div_subset,
                       aes(x = midpoint,
                           y = nuc_diversity)) +
        geom_point(colour = "black",
                   size = 4,
                   alpha = 0.75) +
        xlab("Genomic coordinates") +
        ylab("Nucleotide diversity") +
        geom_vline(xintercept = c(start + 100000,
                                  end - 80000),
                   linetype = "dashed",
                   colour = "blue") +
        geom_hline(yintercept = c(nucdiv_median),
                   linetype = "dashed",
                   colour = "orange") +
        geom_rect(aes(xmin = 764617,
                      xmax = 879060,
                      ymin = -Inf,
                      ymax = Inf),
                  fill = "light blue",
                  alpha = 0.008) +
        theme_bw()
nuc_div_plot <- nuc_div_plot +
        theme(axis.text = element_text(size = 15),
              axis.title = element_text(size = 15,
                                        face = "bold"))
tajima_plot <- ggplot(nuc_div_subset,
                      aes(x = midpoint,
                          y = tajima_d)) +
        geom_point()

## Generate a combined plot:
ggarrange(nsl_plot,
          nuc_div_plot,
          nrow = 2,
          ncol = 1,
          align = "hv",
          labels = c("a",
                     "b"))
## Save:
ggsave(file = paste("results/", output, sep = ""),
       width = 10,
       height = 10)
