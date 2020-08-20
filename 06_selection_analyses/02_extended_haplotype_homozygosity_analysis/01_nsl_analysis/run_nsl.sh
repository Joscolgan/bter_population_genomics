#!/bin/sh

module load r

mkdir input/variant_files
mkdir input/hap
mkdir input/map
mkdir results

## Take input from the command line:
input_vcf=$1

cd input
## Create a file containing list of chromosome/LG names:
grep -v '#' "$input_vcf" | \
cut -f 1 | sort | uniq > chromosome_list.txt

## Using the chromosome list, create individual vcf files for each chromosome:
while read line
do
vcftools --vcf "$input_vcf"  \
--chr "$line" \
--recode \
--out variant_files/"$line"
done < chromosome_list.txt

## The next step would be to create a haplotype file for each individual:
for name in variant_files/*.vcf;
do
abbrev_name="$(echo "$name" | cut -d '.' -f 1,2 - | cut -d '/' -f 2 )";
Rscript ../VCF_to_genomatrix_hap.R "$name" hap/"$abbrev_name".gds hap/"$abbrev_name".hap.out map/"$abbrev_name".map.out;
done

## Generate a genetic map:
## This approach works for nsl where genetic distance between snps is not required - only physical distance:
for name in variant_files/*vcf;
do
new_name="$(echo "$name" | cut -d '.' -f 1,2 | cut -d '/' -f 2 )";
## Run plink
plink --vcf "$name" \
      --double-id \
      --allow-extra-chr \
      --out map/"$new_name";
## Subset columns of interest:
cut -f 1,2,3,4 map/"$new_name".bim > map/"$new_name".tmp;
awk '{ print $1,$2,$4,$4 }' map/"$new_name".tmp | tr ' ' '\t' > map/"$new_name".map;
done

## Remove tmp files:
rm -r vcf/

## Run nsl
for name in hap/*.out;
do
new_name="$(echo "$name" | cut -d '/' -f 2 | cut -d '.' -f 1,2 )";
selscan \
--threads 20 \
--nsl \
--hap "$name" \
--map map/"$new_name".map \
--out ../results/"$new_name";
done

## Normalise output:
norm --nsl --files ../results/*out

cd ../results
for name in *norm;
do
chrom_name="$(echo "$name" | cut -d '.' -f 1 )";
while read line;
do
echo "$chrom_name" >> "$chrom_name".tmp;
done < "$name";
paste "$chrom_name".tmp "$name" >> "$name"_plus_chrom.txt; done 
