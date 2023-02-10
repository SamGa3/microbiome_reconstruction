
#---------------------------------------------------------------------------------------------------------------------------

## TCGA

#---------------------------------------------------------------------------------------------------------------------------

## ALL TISSUES (COAD GBM LUAD LUSC HNSC OV READ SKCM BRCA) 

# ComBat correction corr_plate_id 
# merged_unamb_score_norm

# COAD noFFPE OnlyPrimary NoRiboZeroG
# LUAD noFFPE OnlyPrimary NoRiboZeroG
# LUSC OnlyPrimary NoDupl
# HNSC noFFPE OnlyPrimary
# OV OnlyPrimary NomirVana
# READ OnlyPrimary
# SKCM OnlyPrimary
# BRCA noFFPE OnlyPrimary NoRiboZeroG

# check intra extra of plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction_comparison.Rmd", 
    params=list(
        tissues = c("COAD", "LUAD", "LUSC", "HNSC", "OV", "READ", "SKCM", "BRCA"),
        metadata = list("../../metadata/COAD/COAD_technical_metadata.txt",
                    "../../metadata/LUAD/LUAD_technical_metadata.txt", 
                    "../../metadata/LUSC/LUSC_technical_metadata.txt", 
                    "../../metadata/HNSC/HNSC_technical_metadata.txt", 
                    "../../metadata/OV/OV_technical_metadata.txt", 
                    "../../metadata/READ/READ_technical_metadata.txt", 
                    "../../metadata/SKCM/SKCM_technical_metadata.txt", 
                    "../../metadata/BRCA/BRCA_technical_metadata.txt"),
        taxa_raw = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_corrected = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/HNSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/OV_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/SKCM_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/BRCA_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
                            ),
        batches = rep("plate_id", 8)
    ), 
    output_file = "../../results/technical_batch_effect/COAD_LUAD_LUSC_HNSC_OV_READ_SKCM_BRCA_ComBat_corr_plate_id_bacteria_species_batch_comparison_plate_id.html"
)

rm(list=ls())
gc(full=TRUE)

## COAD and READ 

# ComBat correction corr_plate_id 
# merged_unamb_score_norm

# COAD noFFPE OnlyPrimary NoRiboZeroG
# READ OnlyPrimary

# check intra extra of plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction_comparison.Rmd", 
    params=list(
        tissues = c("COAD", "READ"),
        metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/READ/READ_technical_metadata.txt"),
        taxa_raw = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_corrected = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        batches = rep("read_length", 2)
    ), 
    output_file = "../../results/technical_batch_effect/COAD_READ_ComBat_corr_plate_id_bacteria_species_batch_comparison_read_length.html"
)
