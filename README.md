## Context & Citation:  
This repository contains all scripts used for the population genomics analysis of wild-caught _Bombus terrestris_ males collected across the island of Great Britain.  
This study was developed as part of a NERC-funded project reported implemented in [Yannick Wurm's lab at Queen Mary University of London](https://wurmlab.github.io/).  
Findings are reported in the following manuscript:

Colgan, T.J., Arce, A.N., Gill, R.J., Ramos Rodrigues, A., Kanteh, A., Duncan, E., Li, L., Chittka, L., Wurm, Y.  
__Genomics of adaptation to a changing world in a wild pollinator__, in prep for submission.  

The present repository contains scripts for:  
01_quality_assessment/  
- The quality assessment of raw sequencing data using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).  
02_species_identification/  
- Filtering of raw Illumina FASTQ data (data_analysis/) using the following steps:
- 1) Adaptor identification and removal using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). 
- 2) Interleaving of adaptor removed pairs using [Khmer](https://github.com/dib-lab/khmer). 
- 3) Filtering based on sequence quality using [fastx](http://hannonlab.cshl.edu/fastx_toolkit/).  
- 4) _K_-mer counting and filtering using [Khmer](https://github.com/dib-lab/khmer).  
- 5) Split back into pairs and orphans produced during filtering using [Khmer](https://github.com/dib-lab/khmer).  
- 6) Remove sequences less than 50bp in length using [seqtk](https://github.com/lh3/seqtk).  
- The alignment and variant calling of filtered Illumina data:
- 1) Alignment against the [bumblebee](https://www.ncbi.nlm.nih.gov/assembly/GCF_000214255.1) reference genome was performed using [bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
- The assessment of population structure and admixture (population_structure_analyses/).  
- Investigate signatures of selection, including:  
  - 1) Investigating regions of reduced nucleotide diversity using [PopGenome](https://cran.r-project.org/web/packages/PopGenome/index.html).  
  - 2) Investigating regions of extended haplotype homozygosity using [selscan](https://github.com/szpiech/selscan).  
03_variant_calling/  
04_population_structure_analysis/  
05_selection_analyses/  
06_copy_number_variation_analysis  
- Copy number variation   
07_functional_term_enrichment_analysis/  
- Gene ontology term enrichment analysis using [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html).   
- Functional domain enrichment analysis
- Investigation of genes of interest:  
  - 1) Insecticide response genes  
  - 2) _Wolbachia_ HGT event  
  - 3) Region of low nucleotide diversity on chromosome one  
- Scripts for visualisation of data for manuscript/publication (08_graphics/)  
