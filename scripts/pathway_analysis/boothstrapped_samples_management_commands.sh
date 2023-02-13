#!/bin/bash

# left
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/left_random/COAD_selectedTumor_left_random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/left_random/
    input_name="data/RNAseq/humann_output/left_random/COAD_selectedTumor_left_random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/left_random/COAD_selectedTumor_left_random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# right
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/right_random/COAD_selectedTumor_right_random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/right_random/
    input_name="data/RNAseq/humann_output/right_random/COAD_selectedTumor_right_random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/right_random/COAD_selectedTumor_right_random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# CMS1
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/CMS1_random/COAD_selectedTumor_CMS1_random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/CMS1_random/
    input_name="data/RNAseq/humann_output/CMS1_random/COAD_selectedTumor_CMS1_random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/CMS1_random/COAD_selectedTumor_CMS1_random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# CMS234
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/CMS234_random/COAD_selectedTumor_CMS234_random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/CMS234_random/
    input_name="data/RNAseq/humann_output/CMS234_random/COAD_selectedTumor_CMS234_random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/CMS234_random/COAD_selectedTumor_CMS234_random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# mutation_burden_high
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/mutation_burden_high_random/COAD_selectedTumor_mutation_burden_high_random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/mutation_burden_high_random/
    input_name="data/RNAseq/humann_output/mutation_burden_high_random/COAD_selectedTumor_mutation_burden_high_random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/mutation_burden_high_random/COAD_selectedTumor_mutation_burden_high_random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# mutation_burden_low
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/mutation_burden_low_random/COAD_selectedTumor_mutation_burden_low_random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/mutation_burden_low_random/
    input_name="data/RNAseq/humann_output/mutation_burden_low_random/COAD_selectedTumor_mutation_burden_low_random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/mutation_burden_low_random/COAD_selectedTumor_mutation_burden_low_random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done
