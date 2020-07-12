#!/bin/env bash
##
## Script for reformatting GO term into GeneSet format for use with Gowinda.
##

## Remove all file if exists:
rm gene_list.tmp

## Subset gene ontology terms and reformat:
while read line;
do
  gene_name="$(echo "$line" | cut -f 1 - )";
  go_terms="$(echo "$line" | cut -f 2 - )";
  reformatted="$(echo "$go_terms" | tr ',' '\n' )";
while read term;
do
  echo "$gene_name" >> "$gene_name".tmp;
done < <(echo "$reformatted");
## Combine GO terms and gene names:
paste <(echo "$reformatted") "$gene_name".tmp | sed 's/ /\t/g' >> gene_list.tmp
done < dmel_vs_bter_biomart.input_for_converter.output.txt

## Remove redundant tmp files:
rm LOC*

## For the next step:
## Sort GO terms and create a unique list of terms:
sort -k1,1nr gene_list.tmp | cut -f 1 | uniq > unique_go_term_list.tmp
##
while read unique_go_term;
do
  go_term_of_interest="$(echo "$unique_go_term" | cut -f 1 - )"
  grep "$go_term_of_interest" gene_list.tmp | cut -f 2 - >> "$go_term_of_interest".tmp
  cat "$go_term_of_interest".tmp | tr '\n' ',' >> "$go_term_of_interest".2.tmp
  paste <(echo "$go_term_of_interest") "$go_term_of_interest".2.tmp >> gene_list.2.tmp
  rm "$go_term_of_interest".tmp
  rm "$go_term_of_interest".2.tmp
done < unique_go_term_list.tmp


