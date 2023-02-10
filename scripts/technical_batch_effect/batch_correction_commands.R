
#---------------------------------------------------------------------------------------------------------------------------

## COAD noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

# corr_plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        property = "corr_plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/COAD_selectedTumor_ComBat_batch_correction_corr_plate_id.html"
)

rm(list=ls())
gc(full=TRUE)

# read_length
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        property = "read_length",
        output = "../../data/RNAseq/bacteria/ComBat_read_length/merged_unamb_score_norm/COAD_ComBat_read_length_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/COAD_selectedTumor_ComBat_batch_correction_read_length.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## COAD noFFPE CoupledNormalPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

# corr_plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedCoupledTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        property = "corr_plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedCoupledTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/COAD_selectedCoupledTumorNormal_ComBat_batch_correction_corr_plate_id.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## LUAD noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

# corr_plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/LUAD/LUAD_clinical_metadata.txt", "../../metadata/LUAD/LUAD_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        property = "corr_plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/LUAD_selectedTumor_ComBat_batch_correction_corr_plate_id.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# LUSC OnlyPrimary NoDupl

#---------------------------------------------------------------------------------------------------------------------------

# corr_plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/LUSC/LUSC_clinical_metadata.txt", "../../metadata/LUSC/LUSC_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        property = "corr_plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/LUSC_selectedTumor_ComBat_batch_correction_corr_plate_id.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# HNSC noFFPE OnlyPrimary

#---------------------------------------------------------------------------------------------------------------------------

# corr_plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/HNSC/HNSC_clinical_metadata.txt", "../../metadata/HNSC/HNSC_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        property = "corr_plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/HNSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/HNSC_selectedTumor_ComBat_batch_correction_corr_plate_id.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# OV OnlyPrimary NomirVana

#---------------------------------------------------------------------------------------------------------------------------

# corr_plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/OV/OV_clinical_metadata.txt", "../../metadata/OV/OV_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        property = "corr_plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/OV_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/OV_selectedTumor_ComBat_batch_correction_corr_plate_id.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# READ OnlyPrimary

#---------------------------------------------------------------------------------------------------------------------------

# corr_plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/READ/READ_clinical_metadata.txt", "../../metadata/READ/READ_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        property = "corr_plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/READ_selectedTumor_ComBat_batch_correction_corr_plate_id.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# SKCM OnlyPrimary

#---------------------------------------------------------------------------------------------------------------------------

# corr_plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/SKCM/SKCM_clinical_metadata.txt", "../../metadata/SKCM/SKCM_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        property = "corr_plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/SKCM_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/SKCM_selectedTumor_ComBat_batch_correction_corr_plate_id.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## BRCA noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

# corr_plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/BRCA/BRCA_clinical_metadata.txt", "../../metadata/BRCA/BRCA_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        property = "corr_plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/BRCA_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/BRCA_selectedTumor_ComBat_batch_correction_corr_plate_id.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## IEO

#---------------------------------------------------------------------------------------------------------------------------

# plate_id
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params=list(
        metadata = c("../../metadata/IEO/IEO_clinical_metadata.txt", "../../metadata/IEO/IEO_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        property = "plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/IEO_ComBat_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/IEO_selectedTumor_ComBat_batch_correction_plate_id.html"
)
