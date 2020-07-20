#!/bin/sh




## First extract gene coordinates from GFF:
gff=$1 
nsl_output=$2
output=$3

## Extract gene coordinates from ensembl GFF:
zgrep -v '#' "$gff" | awk '$3=="gene"' | cut -f 1,4,5,9 | cut -d ';' -f 1 | sed 's/ID=gene://' | grep -v 'GL' > input/genes.bed

## Convert chromosome names from ensembl to refseq:
zgrep -v '#' "$gff" | grep -v 'GL' | cut -f 1 | sort | uniq > input/ensembl_chrom_list.txt
zgrep -v '#' "$nsl_output" | cut -f 1 | sort | uniq > input/refseq_chrom_list.txt
paste input/ensembl_chrom_list.txt input/refseq_chrom_list.txt > input/ensembl_to_refseq_chrom_list.txt

while read line
do
ensembl="$(echo "$line" | cut -f 1)"
refseq="$(echo "$line" | cut -f 2)"
echo converting $ensembl to $refseq
sed -i "s/$ensembl/$refseq/g" input/genes.bed
done < input/ensembl_to_refseq_chrom_list.txt

## Create temporary directory:
mkdir -r input/tmp
cd input/tmp

cp -p ../genes.bed .

sed -i 's/\t/-/g' genes.bed

## Create an individual bed file per gene:
while read line
do
echo creating "$line"
gene="$(echo "$line" | cut -d '.' -f 2 | rev | cut -d '-' -f 1 | rev )"
echo "$line" > "$line".bed
sed -i 's/-/\t/g' "$line".bed
~/src/bedtools2/bin/intersectBed -a ../../"$nsl_output" \
-b "$line".bed >> "$line"_overlap.tmp
nsl="$(sort -k4,4gr "$line"_overlap.tmp | head -n 1)"
echo "$gene" "$nsl" >> ../../"$output"
done < genes.bed
