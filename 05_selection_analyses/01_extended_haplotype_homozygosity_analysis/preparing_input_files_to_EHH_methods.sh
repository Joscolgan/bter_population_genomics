module load vcftools
module load plink
module load R

mkdir input
mkdir results

cd input
## Create a symbolic with original VCF: 
ln -s /data/autoScratch/weekly/btw928/2016-10-11_Bombus_population_genomics/results/2017-11-02_terrestris_combined_reanalysis/results/temp/09_final_clean/filtered_temp/02_sorted/2017-11-23-coverage_1/2017-11-26-sample_n41_for_nsl/freebayes_hap0_minQ_1_minaltfrac_0.25_minCov1.Bter_n41.NC_only.minQ20.snps_only.no_hets.sites_low_freq.rare_variants_free.maxMeanDP100.recode.ann.vcf .

## Create a file containing list of chromosome/LG names:
grep -v '#' freebayes_hap0_minQ_1_minaltfrac_0.25_minCov1.Bter_n41.NC_only.minQ20.snps_only.no_hets.sites_low_freq.rare_variants_free.maxMeanDP100.recode.ann.vcf | cut -f 1 | sort | uniq > chromosome_list.txt

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

## Generate a genetic map:
## This approach works for nsl where genetic distance between snps is not required - only physical distance:
for name in vcf/*vcf; do new_name="$(echo "$name" | cut -d '.' -f 1,2 | cut -d '/' -f 2 )";
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
rm -r tmp/
rm -r vcf/

## Run nsl
for name in input/hap/*.out;
do
new_name="$(echo "$name" | cut -d '/' -f 3 | cut -d '.' -f 1,2 )";
singularity exec /data/SBCS-Informatics/singularity/images/ubuntu_trusty-selscan.img selscan \
--nsl \
--hap "$name" \
--map ./input/map/"$new_name".map \
--out results/"$new_name";
done

## For running iHS, a genetic map containg genetic and physical distances between SNPs is required: 
mkdir map_ihs
ln -s ~/autoscratch_monthly_projects/2016-10-11_Bombus_population_genomics/results/2017-11-02_terrestris_combined_reanalysis/results/2017-11-19_run_ldhelmet/*post.txt .

## prepare files for iHS analysis:
for name in *.txt;
do
awk '{ print $1,$3,$1 }' "$name" | tail -n +4 - >> "$name".tmp;
awk '{ print $2,$3,$1 }' "$name" | tail -n 1 >> "$name".tmp;
while read line;
do
new_name="$(echo "$name" | cut -d '.' -f 1,2 )";
echo "$new_name" >> "$new_name".chrom.tmp;
done < "$name".tmp;
paste "$new_name".chrom.tmp "$name".tmp > "$new_name".input_for_genetic_map.txt;
done

## Ensure there are tabs rather than spaces between columns:
for name in *nput_for_genetic_map.txt;
do
sed -i 's/ /\t/g' "$name";
done

## Generate genetic maps for each chromosome:
for name in *input_for_genetic_map.txt;
do
new_name="$(echo "$name" | cut -d '.' -f 1,2 )";
Rscript make_genetic_map.R "$name" "$new_name".ihs_map.txt;
done