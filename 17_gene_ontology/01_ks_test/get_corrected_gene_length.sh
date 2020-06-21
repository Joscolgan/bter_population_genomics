module load seqtk
module load bedtools

## Copy or link file containing gene length to working directory:
cp -p ~/archive-SBCS-Wurmlab/tjcolgan/2017-11-18-Bter_n51_reanalysis/2019-02-ihs_analysis/results/gff_overlap/Bombus_terrestris.Bter_1.0.42.gff3.gene.bed .

## Ensure that the file from ensembl has same as VCF:
cut -f 1 Bombus_terrestris.Bter_1.0.42.gff3.gene.bed | grep -v 'GL' | sort | uniq > chromosome_list.tmp
zgrep -v '#' freebayes_hap0_minQ_1_minaltfrac_0.25_minCov1.Bter_n41.NC_only.minQ20.snps_only.no_hets.sites_low_freq.rare_variants_free.maxMeanDP100.recode.ann.vcf.gz | \
cut -f 1 | uniq > chromosome_list.tmp.2
paste chromosome_list.tmp chromosome_list.tmp.2 > chromosome_list.txt

## Replace:
while read line;
do
old_name="$(echo "$line" | cut -f 1 )"; # Extract original name (i.e. B01)
new_name="$(echo "$line" | cut -f 2 )"; # Replace with name in VCF
sed -i.tmp "s/$old_name/$new_name/g" Bombus_terrestris.Bter_1.0.42.gff3.gene.bed ;
done < chromosome_list.txt 


## For calculating the corrected gene length, three steps will be taken:
## 1) Using individual BED files for each gene, generate a FASTA file for each:
## 2) Using the FASTA files, use seqtk to calculate length and proportion of N bases:
## 3) Subtract the number of N bases from each gene and output.
mkdir tmp
cd tmp
grep -v 'GL' ../Bombus_terrestris.Bter_1.0.42.gff3.gene.bed > Bombus_terrestris.Bter_1.0.42.gff3.gene.bed.tmp
sed -i 's/\t/-/g' Bombus_terrestris.Bter_1.0.42.gff3.gene.bed.tmp

## Generate individual BED files:
while read line
do
echo "$line" > "$line".bed
done < Bombus_terrestris.Bter_1.0.42.gff3.gene.bed.tmp

## Generate FASTA file for each gene:
for name in NC_0157*bed;
do
sed -i 's/-/\t/g' "$name";
fastaFromBed
-fi ../GCF_000214255.1_Bter_1.0_genomic.fna -bed "$name" -fo "$name".fasta; seqtk comp "$name".fasta | \
awk '$14=$2-$9' | tr ' ' '\t' | cut -f 1,14 > "$name".comp;
done

## Combine 'comp' files in the same order as outlined in Bombus_terrestris.Bter_1.0.42.gff3.gene.bed:
while read line;
do
cat "$line".bed.comp >> new_combined.tmp;
done < Bombus_terrestris.Bter_1.0.42.gff3.gene.bed.tmp

## Combine nsl scores and corrected gene lengths:
paste <(sort -k1,1n nsl_scores_per_gene.tmp) <(sort -k1,1n gene_by_corrected_length.txt) | cut -f 1,2,4 > nsl_scores_per_gene.txt
