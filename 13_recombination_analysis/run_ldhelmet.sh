#!/usr/bin/env bash
##############################################################################
# Author: Joe Colgan                   Program: run_LDHelmet.sh
#
# Date: 07/03/2017
#
# Purpose:
# - Parse and create a VCF file for each chromosome
# - 
##############################################################################

## Load modules for analysis
module load vcftools/0.1.15 
module load use.dev ldhelmet

# Read input and from the command line
input=$1  # The first argument ($1) is input
output=$2 # The second argument ($2) is output

##=========================================================
## Step 1: Parse and create a VCF file for each chromosome
##=========================================================
echo 'Step 1: Parse and create a VCF file for each chromosome'

## Taking the input VCF from the command line
## Change the status of unphased to phased:
sed -i 's/\//\|/g' "$input"

## Parse chromosome information from the first column to allow
## for subsetting:
grep -v '#' "$input" | cut -f 1 - | sort | uniq > chromosome_list.txt

## Using the chromosome_list.txt, generate an individual VCF file 
## for each chromosome:
while read line
do
  vcftools --vcf "$input" \
           --chr "$line" \
           --recode \
           --out ./results/"$line"
done < chromosome_list.txt

echo 'Step 1: Parse and create a VCF file for each chromosome: Complete'

##=========================================================
## Step 2: Generate input file for LDHelmet
##=========================================================
echo 'Step 2: Generating input files for LDHelmet'

for name in results/*.recode.vcf
do
new_name="$(echo "$name" | cut -d '/' -f 2 | cut -d '.' -f 1,2 - )"
  vcftools --vcf "$name" \
           --ldhelmet \
           --chr "$new_name" \
           --out "$name"
done     

echo 'Step 2: Generating input files for LDHelmet - complete'

##
for name in results/*snps; 
do 
grep '>' "$name" >> results/hap_names.txt; 
done

##
sed -i 's/>//g' results/hap_names.txt

##
grep -v '\-1' results/hap_names.txt | sort | uniq - > results/hap_names.uniq.txt

##
for name in results/*.snps; 
do 
~/src/seqtk/seqtk subseq "$name" results/hap_names.uniq.txt > "$name".hap; 
done

##=========================================================
## Step 3: Running LDHelmet
##=========================================================
echo 'Step 3: Running LDHelmet'

for name in results/*.snps.hap
do
ldhelmet find_confs \
          --num_threads 10 \
          -w 50 \
          -o "$name".conf \
          "$name"
done

echo 'Step 4: Generating genetic table'
for name in results/*.conf
do
ldhelmet table_gen \
         --num_threads 10 \
          -c "$name" \
          -t 0.00075 \
          -r 0.0 0.1 10.0 1.0 100.0 \
          -o "$name".lk
done

echo 'Step 4: Generating genetic table - complete'

echo 'Step 5: Generating Pade coefficients'

for name in results/*.conf
do
ldhelmet pade \
         --num_threads 10 \
          -c "$name" \
          -t 0.00075 \
          -x 11 \
          -o "$name".pade
done

echo 'Step 5: Generating Pade coefficients - complete'

echo 'Step 6: Generating population-scaled recombination rates'

for name in results/*conf
do
new_name="$( echo "$name" | cut -d '.' -f 1,2,3,4,5 - )"
ldhelmet rjmcmc \
        --num_threads 10 \
        -w 50 \
        -l "$name".lk \
        -p "$name".pade \
        --snps_file "$new_name".snps.hap \
        --pos_file  "$new_name".pos \
        -b 5.0 \
        --burn_in 100000 \
        -n 1000000 \
        -o "$new_name".ldhelmet.post
done

echo 'Step 6: Generating population-scaled recombination rates - complete'

echo 'Step 7: Converting from post to text'

for name in results/*post
do
ldhelmet post_to_text \
  -m -p 0.025 -p 0.50 -p 0.975 \
  -o  "$name".txt \
  "$name"
done

echo 'Step 7: Converting from post to text - complete'
