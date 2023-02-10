#!/bin/bash

# left
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/random_subset_left/COAD_selectedTumor_left_55random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/random_subset_left/
    input_name="data/RNAseq/humann_output/random_subset_left/COAD_selectedTumor_left_55random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/random_subset_left/COAD_selectedTumor_left_55random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# right
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/random_subset_right/COAD_selectedTumor_right_72random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/random_subset_right/
    input_name="data/RNAseq/humann_output/random_subset_right/COAD_selectedTumor_right_72random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/random_subset_right/COAD_selectedTumor_right_72random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# CMS1
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/random_subset_CMS1/COAD_selectedTumor_CMS1_17random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/random_subset_CMS1/
    input_name="data/RNAseq/humann_output/random_subset_CMS1/COAD_selectedTumor_CMS1_17random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/random_subset_CMS1/COAD_selectedTumor_CMS1_17random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# CMS234
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/random_subset_CMS234/COAD_selectedTumor_CMS234_86random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/random_subset_CMS234/
    input_name="data/RNAseq/humann_output/random_subset_CMS234/COAD_selectedTumor_CMS234_86random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/random_subset_CMS234/COAD_selectedTumor_CMS234_86random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# mutation_burden_high
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/random_subset_mutation_burden_high/COAD_selectedTumor_mutation_burden_high_17random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/random_subset_mutation_burden_high/
    input_name="data/RNAseq/humann_output/random_subset_mutation_burden_high/COAD_selectedTumor_mutation_burden_high_17random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/random_subset_mutation_burden_high/COAD_selectedTumor_mutation_burden_high_17random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# mutation_burden_low
for i in {1..50..1}
do
    input_name="data/RNAseq/humann_output/random_subset_mutation_burden_low/COAD_selectedTumor_mutation_burden_low_95random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/random_subset_mutation_burden_low/
    input_name="data/RNAseq/humann_output/random_subset_mutation_burden_low/COAD_selectedTumor_mutation_burden_low_95random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/random_subset_mutation_burden_low/COAD_selectedTumor_mutation_burden_low_95random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done
