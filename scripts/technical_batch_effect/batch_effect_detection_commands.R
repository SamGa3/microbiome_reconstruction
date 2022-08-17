
#---------------------------------------------------------------------------------------------------------------------------

## TCGA

#---------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
gc(full=TRUE)

## ALL TISSUES (COAD GBM LUAD LUSC HNSC OV READ) 

# merged_unamb_score_norm

# COAD noFFPE OnlyPrimary NoRiboZeroG
# LUAD noFFPE OnlyPrimary NoRiboZeroG
# LUSC OnlyPrimary NoDupl
# HNSC noFFPE OnlyPrimary
# OV OnlyPrimary NomirVana
# READ OnlyPrimary
# SKCM OnlyPrimary

rmarkdown::render("scripts/technical_batch_effect/batch_detection.Rmd", 
  params=list(
    tissues = c("COAD", "GBM", "LUAD", "LUSC", "HNSC", "OV", "READ", "SKCM"),
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/GBM/GBM_technical_metadata.txt",
                  "../../metadata/LUAD/LUAD_technical_metadata.txt", "../../metadata/LUSC/LUSC_technical_metadata.txt", 
                  "../../metadata/HNSC/HNSC_technical_metadata.txt", "../../metadata/OV/OV_technical_metadata.txt", 
                  "../../metadata/READ/READ_technical_metadata.txt", "../../metadata/SKCM/SKCM_technical_metadata.txt"),
    new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
    taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
                )
  ), 
  output_file = "../../results/technical_batch_effect/COAD_LUAD_LUSC_HNSC_OV_READ_SKCM_bacteria_species_merged_unamb_score_norm_batch_detection.html"
)

rm(list=ls())
gc(full=TRUE)

## COAD 

# merged_unamb_score_norm merged_unamb_norm merged_score_norm merged_ambig_norm

# COAD noFFPE OnlyPrimary NoRiboZeroG

rmarkdown::render("scripts/technical_batch_effect/batch_detection.Rmd", 
  params=list(
    tissues = c("merged_unamb_score_norm", "merged_unamb_norm", "merged_ambig_norm"),
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"),
    new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
    taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
              "../../data/RNAseq/bacteria/raw/merged_unamb_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_norm.txt",
              "../../data/RNAseq/bacteria/raw/merged_ambig_norm/COAD/COAD_selectedTumor_bacteria_species_merged_ambig_norm.txt")
  ), 
  output_file = "../../results/technical_batch_effect/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm_merged_unamb_norm_merged_ambig_norm_batch_detection.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## IEO

#---------------------------------------------------------------------------------------------------------------------------

# merged_unamb_score_norm

rmarkdown::render("scripts/technical_batch_effect/batch_detection.Rmd", 
  params=list(
    tissues = c("IEO"),
    metadata = list(c("../../metadata/IEO/IEO_technical_metadata.txt")),
    join = list("rows"),
    exceptions = c("plate_id"),
    taxa = list("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumor_bacteria_species_merged_unamb_score_norm.txt")
    ), 
  output_file = "../../results/technical_batch_effect/IEO_selectedTumor_bacteria_species_merged_unamb_score_norm_batch_detection.html"
)

rm(list=ls())
gc(full=TRUE)
