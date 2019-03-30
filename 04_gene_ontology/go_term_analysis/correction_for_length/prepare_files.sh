mkdir gene_coordinates
mkdir enriched_sites
cd gene_coordinates/
ln -s ../../Bombus_terrestris.Bter_1.0.42.gff3.gene.bed .
cd ../
cd enriched_sites/
cut -f 1,2,3,7 nsl_n41_unfiltered.txt > nsl_n41_unfiltered_corrected_pvalues.bed
cut -f 1,2,3,6 nsl_n41_unfiltered.txt > nsl_n41_unfiltered_raw_pvalues.bed
cd ../
mkdir tmp_raw_values
mkdir tmp_corrected_values

cd tmp_raw_values
mkdir genes
mkdir sites
cd genes
ln -s ../../gene_coordinates/Bombus_terrestris.Bter_1.0.42.gff3.gene.bed .
cut -f 1  Bombus_terrestris.Bter_1.0.42.gff3.gene.bed | uniq | grep 'B' > old_chromosome_names.tmp
cut -f 1 ../../enriched_sites/nsl_n41_unfiltered.txt | uniq | paste old_chromosome_names.tmp - > new_chromosome_names.txt
## For testing:
cat Bombus_terrestris.Bter_1.0.42.gff3.gene.bed > test.bed

while read line;
do
old_chrom="$(echo "$line" | cut -f 1 )";
new_chrom="$(echo "$line" | cut -f 2 )";
echo "$old_chrom";
echo "$new_chrom";
sed -i "s/$old_chrom/$new_chrom/g" test.bed ;
done < new_chromosome_names.txt 

mkdir tmp
cd tmp
cp -p ../test.bed .
sed -i 's/\t/-/g' test.bed
while read line; do echo "$line" > "$line".bed; done < test.bed 
for name in *.bed; do sed -i 's/-/\t/g' "$name"; done
for name in *.bed; do intersectBed -a ../../../enriched_sites/nsl_n41_unfiltered_corrected_pvalues.bed -b "$name" > "$name"_nsl_site_overlap.tmp; done
for name in *.bed; do intersectBed -a ../../../enriched_sites/nsl_n41_unfiltered_raw_pvalues.bed -b "$name" > "$name"_nsl_site_overlap_raw.tmp; done
