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

##############################################################################
## Step One: Convert input VCF to plink file format
##############################################################################

## Print to console:
echo "Step One: Converting VCF to plink file format"

## Use vcftools to convert VCF to PED file format:
vcftools --vcf "$input" \
         --plink \
         --out "$input".plink

## Use plink to convert PED file to BED files:
plink --file "$input".plink \
         --maf 0.1 \
         --make-bed \
         --out "$input".plink

## Print to console:
echo "Step One - complete!"

##############################################################################
## Step Two: Convert input VCF to plink file format
##############################################################################

## Print to console:
echo "Step Two: Running admixture"

## Run admixture with 100 interactions for 10 potential populations
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15;
do
../src/admixture_linux-1.3.0/admixture --cv "$input".plink.bed -j20  $K | tee log${K}.out;
done

## Print to console:
echo "Step Two - complete!"
