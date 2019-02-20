#!/bin/env bash
## Load modules
module load vcftools
module load plink
module load R

## Create directories for analysis:
mkdir input
mkdir results
mkdir code

cd input
mkdir vcf
mkdir hap
mkdir map
## Create a symbolic with original VCF:
ln -s /data/autoScratch/weekly/btw928/2016-10-11_Bombus_population_genomics/results/2017-11-02_terrestris_combined_reanalysis/results/temp/09_final_clean/filtered_temp/02_sorted/2017-11-23-coverage_1/201$

## Create a file containing list of chromosome/LG names:
grep -v '#' freebayes_hap0_minQ_1_minaltfrac_0.25_minCov1.Bter_n41.NC_only.minQ20.snps_only.no_hets.sites_low_freq.rare_variants_free.maxMeanDP100.recode.ann.vcf | cut -f 1 | sort | uniq > chromosome_lis$

## Using the chromosome list, create individual vcf files for each chromosome:
while read line
do
vcftools --vcf freebayes_hap0_minQ_1_minaltfrac_0.25_minCov1.Bter_n41.NC_only.minQ20.snps_only.no_hets.sites_low_freq.rare_variants_free.maxMeanDP100.recode.ann.vcf \
--chr "$line" \
--recode \
--out vcf/"$line"
done < chromosome_list.txt

## The next step would be to create a haplotype file for each individual:
for name in vcf/*.vcf;
do
abbrev_name="$(echo "$name" | cut -d '.' -f 1,2 - | cut -d '/' -f 2 )";
Rscript VCF_to_genomatrix_SNPs.R "$name" ./hap/"$abbrev_name".gds ./hap/"$abbrev_name".hap.out;
done
