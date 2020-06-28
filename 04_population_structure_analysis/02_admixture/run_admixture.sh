#!/usr/bin/env bash
##############################################################################
# Author: Joe Colgan                   Program: run_admixture.sh
#
# Date: 24/08/2017
#
# Purpose:
# The script takes a VCF file as input.
# The first step involves the conversion of VCF to PLINK file format.
# The second step runs ADMIXTURE for a number of used defined populations.
#
##############################################################################

## Load modules:
module load vcftools
module load plink

## Take input from command line:
input=$1 ## The first argument on the command line will be taken as input

abbrev_name="$(echo "$input" | cut -d '.' -f 1,2 )"

## Ensure output directory is created:
mkdir -p results/admixture/"$abbrev_name"

##############################################################################
## Step One: Convert input VCF to plink file format
##############################################################################

## Print to console:
echo "Step One: Converting VCF to plink file format"

## Use vcftools to convert VCF to PED file format:
vcftools --vcf results/"$input" \
         --plink \
         --out results/admixture/"$abbrev_name"/"$input".plink

## Use plink to convert PED file to BED files:
plink --file results/admixture/"$abbrev_name"/"$input".plink \
         --maf 0.1 \
         --make-bed \
         --out results/admixture/"$abbrev_name"/"$input".plink

## Print to console:
echo "Step One - complete!"

##############################################################################
## Step Two: Convert input VCF to plink file format
##############################################################################

## Print to console:
echo "Step Two: Running admixture"

## Run admixture with 100 interactions for 10 potential populations
for K in {1..20};
do
admixture --cv results/admixture/"$abbrev_name"/"$input".plink.bed \
-j20 $K | tee results/admixture/"$abbrev_name"/log${K}.out;
done

## Print to console:
echo "Step Two - complete!"
