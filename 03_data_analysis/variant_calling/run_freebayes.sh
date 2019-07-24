freebayes \
-f ./data/refseq_files/GCF_000214255.1_Bter_1.0_genomic.fna \
--bam-list bam_list.txt \
--ploidy 2 \
--report-genotype-likelihood-max \
--use-mapping-quality \
--genotype-qualities \
--use-best-n-alleles 4 \
--haplotype-length 0 \
--min-base-quality 3 \
--min-mapping-quality 1 \
--min-alternate-fraction 0.25 \
--min-coverage 1 \
--use-reference-allele > results/freebayes_hap0_minQ_1_minaltfrac_0.25_minCov1.vcf 

