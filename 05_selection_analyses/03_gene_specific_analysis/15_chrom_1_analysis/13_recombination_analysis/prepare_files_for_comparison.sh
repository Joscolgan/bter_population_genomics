mkdir input
ln -s ../../2019-08-09b_ldhelmet_run_n46/results/*post.txt .
mkdir results


for name in input/*tmp;
do
chrom_name="$(echo "$name" | cut -d '/' -f 2 | cut -d '.' -f 1,2 )";
while read line; do echo "$chrom_name" >> input/"$chrom_name".tmp;
done < "$name";
done

## Create a temporary file with the name of the chromosome:
for name in input/*tmp;
do
chrom_name="$(echo "$name" | cut -d '/' -f 2 | cut -d '.' -f 1,2 )";
echo "$chrom_name";
paste input/"$chrom_name".tmp "$name" > input/"$chrom_name"_combined.txt;
done

## Generate a combined file:
cat input/*combined.txt > input/genome_wide_pop_scaled_recombination_rates.txt
~/src/bedtools2/bin/intersectBed -a ./input/genome_wide_pop_scaled_recombination_rates.txt -b ./data/region_of_interest.bed > ./input/region_of_interest.txt
