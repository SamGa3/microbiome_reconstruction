#!/bin/bash

# left
for i in {1..30..1}
do
    input_name="data/RNAseq/humann_output/random_subset_left/COAD_selectedTumor_left_50random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/random_subset_left/
    input_name="data/RNAseq/humann_output/random_subset_left/COAD_selectedTumor_left_50random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/random_subset_left/COAD_selectedTumor_left_50random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done

# right
for i in {1..30..1}
do
    input_name="data/RNAseq/humann_output/random_subset_right/COAD_selectedTumor_right_50random_${i}_pathabundance.tsv" 
    humann_split_stratified_table --input $input_name --output data/RNAseq/humann_output/random_subset_right/
    input_name="data/RNAseq/humann_output/random_subset_right/COAD_selectedTumor_right_50random_${i}_pathabundance_unstratified.tsv" 
    output_name="data/RNAseq/humann_output/random_subset_right/COAD_selectedTumor_right_50random_${i}_pathabundance_unstratified_cpm.tsv" 
    humann_renorm_table --input $input_name \
                    --output $output_name \
                    --units cpm \
                    --update-snames
done


