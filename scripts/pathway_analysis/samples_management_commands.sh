#!/bin/bash

# left
humann_split_stratified_table --input data/RNAseq/humann_output/left/COAD_selectedTumor_left_pathabundance.tsv --output data/RNAseq/humann_output/left/
humann_renorm_table --input data/RNAseq/humann_output/left/COAD_selectedTumor_left_pathabundance_unstratified.tsv \
                    --output data/RNAseq/humann_output/left/COAD_selectedTumor_left_pathabundance_unstratified_cpm.tsv \
                    --units cpm \
                    --update-snames
# right
humann_split_stratified_table --input data/RNAseq/humann_output/right/COAD_selectedTumor_right_pathabundance.tsv --output data/RNAseq/humann_output/right/
humann_renorm_table --input data/RNAseq/humann_output/right/COAD_selectedTumor_right_pathabundance_unstratified.tsv \
                    --output data/RNAseq/humann_output/right/COAD_selectedTumor_right_pathabundance_unstratified_cpm.tsv \
                    --units cpm \
                    --update-snames
                    
# CMS1
humann_split_stratified_table --input data/RNAseq/humann_output/CMS1/COAD_selectedTumor_CMS1_pathabundance.tsv --output data/RNAseq/humann_output/CMS1/
humann_renorm_table --input data/RNAseq/humann_output/CMS1/COAD_selectedTumor_CMS1_pathabundance_unstratified.tsv \
                    --output data/RNAseq/humann_output/CMS1/COAD_selectedTumor_CMS1_pathabundance_unstratified_cpm.tsv \
                    --units cpm \
                    --update-snames
# CMS234
humann_split_stratified_table --input data/RNAseq/humann_output/CMS234/COAD_selectedTumor_CMS234_pathabundance.tsv --output data/RNAseq/humann_output/CMS234/
humann_renorm_table --input data/RNAseq/humann_output/CMS234/COAD_selectedTumor_CMS234_pathabundance_unstratified.tsv \
                    --output data/RNAseq/humann_output/CMS234/COAD_selectedTumor_CMS234_pathabundance_unstratified_cpm.tsv \
                    --units cpm \
                    --update-snames

# mutation_burden_high
humann_split_stratified_table --input data/RNAseq/humann_output/mutation_burden_high/COAD_selectedTumor_mutation_burden_high_pathabundance.tsv --output data/RNAseq/humann_output/mutation_burden_high/
humann_renorm_table --input data/RNAseq/humann_output/mutation_burden_high/COAD_selectedTumor_mutation_burden_high_pathabundance_unstratified.tsv \
                    --output data/RNAseq/humann_output/mutation_burden_high/COAD_selectedTumor_mutation_burden_high_pathabundance_unstratified_cpm.tsv \
                    --units cpm \
                    --update-snames
# mutation_burden_low
humann_split_stratified_table --input data/RNAseq/humann_output/mutation_burden_low/COAD_selectedTumor_mutation_burden_low_pathabundance.tsv --output data/RNAseq/humann_output/mutation_burden_low/
humann_renorm_table --input data/RNAseq/humann_output/mutation_burden_low/COAD_selectedTumor_mutation_burden_low_pathabundance_unstratified.tsv \
                    --output data/RNAseq/humann_output/mutation_burden_low/COAD_selectedTumor_mutation_burden_low_pathabundance_unstratified_cpm.tsv \
                    --units cpm \
                    --update-snames
