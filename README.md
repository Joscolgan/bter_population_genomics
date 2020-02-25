## Context & Citation:  
This repository contains scripts used for the population genomics analyses outlined the manuscript:   
__Signatures of ongoing selection in response to environmental pressures in a wild population of the buff-tailed bumblebee__
Colgan TJ, Arce AN, Gill RJ, Ramos Rodrigues A, Kanteh A, Li L, Chittka L, Wurm Y.  

The present repository contains scripts for:  
- The quality assessment of raw sequencing data using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).  
- Filtering of raw Illumina FASTQ data (data_analysis/) using the following steps:
- 1) Adaptor identification and removal using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). 
- 2) Interleaving of adaptor removed pairs using [Khmer](https://github.com/dib-lab/khmer). 
- 3) Filtering based on sequence quality using [fastx](http://hannonlab.cshl.edu/fastx_toolkit/).  
- 4) K-mer counting and filtering using [Khmer](https://github.com/dib-lab/khmer).  
- 5) Split back into pairs and orphans produced during filtering using [Khmer](https://github.com/dib-lab/khmer).  
- 6) Remove sequences less than 50bp in length using [seqtk](https://github.com/lh3/seqtk).  
- The alignment and variant calling of filtered Illumina FASTQ data (data_analysis/)  
- The assessment of population structure and admixture (population_structure_analyses/)  
- Investigate signatures of selection, including:  
  - Reduced nucleotide diversity  (selection_analysis/reduced_nucleotide_diversity)
  - Extended haplotype homozygosity  (selection_analysis/extended_haplotype_homozygosity)
  - Copy number variation (selection_analysis/copy_number_variation)
- Phylogenetic analysis of candidate genes:  
- Functional domain and gene ontology term enrichment (functional_term_enrichment/)  
- Scripts for visualisation of data for manuscript/publication (manuscript_plots/)  
