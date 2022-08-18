#!/bin/bash

cd /microbiome_reconstruction

# change mod
chmod u+wrx scripts/pathway_analysis/boothstrapped_samples_management_commands.sh

# Untar score tables
tar -xvf data/RNAseq/bacteria/raw/score/COAD_bacteria_species_score.tar.xz --directory data/RNAseq/bacteria/raw/score
tar -xvf data/RNAseq/bacteria/raw/score/OV_bacteria_species_score.tar.xz --directory data/RNAseq/bacteria/raw/score
tar -xvf data/RNAseq/bacteria/raw/score/LUAD_bacteria_species_score.tar.xz --directory data/RNAseq/bacteria/raw/score
tar -xvf data/RNAseq/bacteria/raw/score/HNSC_bacteria_species_score.tar.xz --directory data/RNAseq/bacteria/raw/score
tar -xvf data/RNAseq/bacteria/raw/score/LUSC_bacteria_species_score.tar.xz --directory data/RNAseq/bacteria/raw/score
tar -xvf data/RNAseq/bacteria/raw/score/SKCM_bacteria_species_score.tar.xz --directory data/RNAseq/bacteria/raw/score

### Create microbiome_reconstruction structure of folders

# data folder
mkdir -p data/pathseq_tools
mkdir -p data/RNAseq/pathseq_output/ERR2756905_out
mkdir -p data/RNAseq/pathseq_output/ERR2756906_out
mkdir -p data/RNAseq/pathseq_output/ERR2756907_out
mkdir -p data/RNAseq/humann_output/random_subset_left
mkdir -p data/RNAseq/humann_output/random_subset_right
mkdir -p data/RNAseq/humann_output/example/ERR2756905_out
mkdir -p data/RNAseq/humann_output/example/ERR2756906_out
mkdir -p data/RNAseq/humann_output/example/ERR2756907_out
mkdir -p data/RNAseq/bacteria/ComBat_read_length/merged_unamb_score_norm
mkdir -p data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD
mkdir -p data/RNAseq/bacteria/raw/merged_ambig_norm/COAD
mkdir -p data/RNAseq/bacteria/raw/merged_score_norm/COAD
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_norm/COAD
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO
mkdir -p data/RNAseq/bacteria/raw/merged_unamb_score_norm/example

# result folder
mkdir -p results/comparison
mkdir -p results/filters
mkdir -p results/identification_related_species/tables
mkdir -p results/ml/images
mkdir -p results/pathway_analysis/bootstrapped_samples
mkdir -p results/property_association/bacteria_species/raw/merged_unamb_score_norm/all/images
mkdir -p results/property_association/bacteria_species/raw/merged_unamb_score_norm/COAD
mkdir -p results/property_association/bacteria_species/raw/merged_unamb_score_norm/GBM/tables
mkdir -p results/property_association/bacteria_species/ComBat_batch_corrected_read_length/merged_unamb_score_norm/COAD
mkdir -p results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables
mkdir -p results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/GBM/tables
mkdir -p results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/tables
mkdir -p results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/tables
mkdir -p results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/tables
mkdir -p results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/tables
mkdir -p results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/tables
mkdir -p results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/tables
mkdir -p results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/IEO
mkdir -p results/survival_analysis
mkdir -p results/technical_batch_effect
