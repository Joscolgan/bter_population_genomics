#!/env/usr/bin Rscript
####################################################################
##
## Name: Joe Colgan               Program: manhattan_plottr.R
##
## Date: 02-06-2019
##
## Purpose:
## This script takes two tab-delimited files as input.
## 1) The first contains five columns:
## - Chromosome name (CHR)
## - SNP position (BP)
## - nsl score (P)
## - variant annotation (variant)
## - gene name (locus)
## 2) the second contains four:
## - CHR
## - BP
## - P
## - locus
## The script takes these inputs and generates a manhattan plot, 
## highlighting the SNPs provided in the second input file. It outputs
## a manhattan plot, which is saved as a .png file but can be
## easily changed to a .jpeg or .pdf.
##
###################################################################

## Load libraries:
libraries <- c("devtools",
               "ggplot2",
               "PopGenome",
               "qqman", 
               "reshape",
               "ggpubr",
               "zoo",
               "stringr",
               "dplyr",
               "ggrepel")
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
input_variants <- args[1]                         # Assign first argument as input
input_peaks <- args[2]                        # Assign second argument as input
output <- args[3]
message(paste("will output to", output))

## Print to console:
print("Step one - loading variant data...")
## Read in variants:
variants <- read.table(file = input_variants,
                       col.names = c("CHR",
                                     "BP",
                                     "P",
                                     "variant",
                                     "locus"))

variants$SNP<- rownames(variants)
variants$P <- abs(variants$P)

## Convert CHR to factor and sort levels by correct order:
variants$CHR <- factor(variants$CHR, levels = unique(variants$CHR))

## Print to console:
print("Step one - variant data loaded!")
## Print to console:
print("Step two - loading peak information...")

## Extract top 50 peaks of interest:
peaks <- read.table(file = input_peaks,
                    col.names = c("CHR",
                                  "BP",
                                  "P",
                                  "locus"))

## Print to console:
print("Step two - peak information loaded!")

peaks_tmp <- subset(variants,
                    BP %in% peaks$BP)

## Create SNP list of interest:
snps_of_interest <- peaks_tmp$SNP

variants$CHR <- gsub("NC_015762.1", 1, variants$CHR)
variants$CHR <- gsub("NC_015763.1", 2, variants$CHR)
variants$CHR <- gsub("NC_015764.1", 3, variants$CHR)
variants$CHR <- gsub("NC_015765.1", 4, variants$CHR)
variants$CHR <- gsub("NC_015766.1", 5, variants$CHR)
variants$CHR <- gsub("NC_015767.1", 6, variants$CHR)
variants$CHR <- gsub("NC_015768.1", 7, variants$CHR)
variants$CHR <- gsub("NC_015769.1", 8, variants$CHR)
variants$CHR <- gsub("NC_015770.1", 9, variants$CHR)
variants$CHR <- gsub("NC_015771.1", 10, variants$CHR)
variants$CHR <- gsub("NC_015772.1", 11, variants$CHR)
variants$CHR <- gsub("NC_015773.1", 12, variants$CHR)
variants$CHR <- gsub("NC_015774.1", 13, variants$CHR)
variants$CHR <- gsub("NC_015775.1", 14, variants$CHR)
variants$CHR <- gsub("NC_015776.1", 15, variants$CHR)
variants$CHR <- gsub("NC_015777.1", 16, variants$CHR)
variants$CHR <- gsub("NC_015778.1", 17, variants$CHR)
variants$CHR <- gsub("NC_015779.1", 18, variants$CHR)

## Print to console:
print("Step three - formatting data for plots...")

don <- variants %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(variants, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%

  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% snps_of_interest, "yes", "no")) %>%
  mutate( is_annotate=ifelse(-log10(P)>-0.40824, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

## Print to console:
print("Step three - data formatted!")
## Print to console:
print("Step four - generating manhattan plot...")

# Make the plot
ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR,
                        breaks = axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis

    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +

    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

## Print to console:
print("Step four - manhattan plot generated!")
## Print to console:
print("Step five - saving manhattan plot to file...")

## Save the plot to file:
ggsave(file = output,
       height = 10,
       width = 10)

## Print to console:
print("Step five - plot saved!")

