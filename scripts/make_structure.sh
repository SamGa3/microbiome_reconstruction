#!/bin/bash

cd /microbiome_reconstruction

# Untar score tables
tar -xvzf data/RNAseq/bacteria/raw/score/COAD_bacteria_species_score.tar.xz
tar -xvzf data/RNAseq/bacteria/raw/score/OV_bacteria_species_score.tar.xz
tar -xvzf data/RNAseq/bacteria/raw/score/LUAD_bacteria_species_score.tar.xz
tar -xvzf data/RNAseq/bacteria/raw/score/HNSC_bacteria_species_score.tar.xz
tar -xvzf data/RNAseq/bacteria/raw/score/LUSC_bacteria_species_score.tar.xz
tar -xvzf data/RNAseq/bacteria/raw/score/SKCM_bacteria_species_score.tar.xz

### Create microbiome_reconstruction structure of folders

# data folder
mkdir -p data/pathseq_tools
mkdir -p data/humann_tools/chocophlan
mkdir -p data/humann_tools/uniref
mkdir -p data/RNAseq/humann_output/random_subset_left
mkdir -p data/RNAseq/humann_output/random_subset_right
mkdir -p data/RNAseq/humann_output/example
mkdir -p data/RNAseq/ComBat_read_length/merged_unamb_score_norm
mkdir -p data/RNAseq/ComBat_plate_id/merged_unamb_score_norm

# result folder
mkdir -p results/comparison
mkdir -p results/filters
mkdir -p results/identification_related_species/tables
mkdir -p results/ml
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
