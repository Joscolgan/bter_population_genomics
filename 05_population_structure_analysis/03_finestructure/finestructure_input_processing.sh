#!/usr/bin/env bash  
## Processing of data for running fineSTRUCTURE:
module load plink

## Run plink: Plink converts VCF into a ped and map file. 
## The output can be used by chromopainter to generate formatted data for input into fs.

for name in NC*vcf
do
plink --vcf "$name" \
      --double-id \
      --recode12 \
      --allow-extra-chr \
      --out "$name"
done

## Run plink2chrompainter.pl:
for name in NC*vcf
do
~/src/fs-2.0.7/fs-2.0.7/scripts/plink2chromopainter.pl \
-p="$name".ped \
-m="$name".map \
-o="$name".phase_file
done

## The genetic map information was obtained from the output of LDHelmet:
for name in *post.txt; do tail -n +4 "$name" > "$name".tmp; done

## Make the recombination map
for name in *.tmp; do awk '{ print $1,$3 }' "$name" > "$name".genetic_map.txt; 
tail -n 1 "$name" | awk '{ print $2, $3 }' - >> "$name".genetic_map.txt;
done

## This results in the generation of information for the phasefile:
for name in *genetic_map.txt
do
Rscript make_genetic_map.R "$name" "$name".recomb.txt
sed -i 's/"//g' "$name".recomb.txt
done 

## Run fineSTRUCTURE:
for name in *.recomb.txt
do
new_name="$(echo "$name | cut -d '.' -f 1,2 - )"
fs "$new_name" -idfile idfile.txt -phasefiles "$new_name".recode.vcf.phase_file \
-recombfiles "$name" -go
done

## Run for all chromosomes together:
~/src/fs-2.0.7/fs bter_fs_analysis -idfile ../data_ids.txt --phasefiles *.phase_file -recombfiles *.genetic_map.txt -go
