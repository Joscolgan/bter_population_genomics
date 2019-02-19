#!/bin/env bash
java -Xmx4G -jar src/Gowinda-1.12.jar \
--output-file ./results/test_run.txt \
--snp-file ./input/total_snp_file.txt \
--candidate-snp-file ./input/snps_of_interest.txt \
--gene-set-file ./input/mart_export_dmel_bter_go_terms_converted.new.txt \
--annotation-file ./data/ensembl/Bombus_terrestris.Bter_1.0.42.new.gtf \
--simulations 10000000 \
--min-significance 1 \
--min-genes 10 \
--gene-definition gene \
--threads 20
