#### TJ Colgan 2019

This repository contains scripts for preparing input files and running extended haplotype homozygosity tests as implemented by [selscan](https://github.com/szpiech/selscan).  
Specifically, the scripts calculate:  
1. The integrated haplotype score (iHS) score.  
2. nSL score.  

### Preparing input files:  
For both tests, a haplotype file and map file are required.  
The haplotype file consists of a genotype matrix whereby each column represents an individual variant and row represents an individual sample.  
The haplotype file can be generated for each individual chromosome using ```VCF_to_geno```.  
The input for the script is a VCF.  
