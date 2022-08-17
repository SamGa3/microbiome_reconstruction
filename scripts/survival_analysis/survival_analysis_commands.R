### COMMANDS

#---------------------------------------------------------------------------------------------------------------------------

## COX

#---------------------------------------------------------------------------------------------------------------------------

## COAD noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
gc(full=TRUE)

# PC1 2 3 4 5 6
rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
    params=list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt",
                        "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_pca_tab.txt"),
        join = c("columns"),
        taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        surv_data = "../../metadata/COAD/COAD_cBioPortal_overall_survival_firehose.txt",
        survival_analysis = c("OS_YEARS", "patient_status", "OS"),
        categorical_covariates = "",
        values_not_considered = "",
        timerange_cat = "",
        numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
        timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5))
    ), 
    output_file = "../../results/survival_analysis/COAD_selectedTumor_COX_OS.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
    params=list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt",
                        "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_pca_tab.txt"),
        join = c("columns"),
        taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        surv_data = "../../metadata/COAD/COAD_cBioPortal_disease_free_survival_firehose.txt",
        survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
        categorical_covariates = "",
        values_not_considered = "",
        timerange_cat = "",
        numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
        timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5))
    ), 
    output_file = "../../results/survival_analysis/COAD_selectedTumor_COX_DFS.html"
)

rm(list=ls())
gc(full=TRUE)

# PC4

rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
    params=list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt",
                        "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_pca_tab.txt"),
        join = c("columns"),
        taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        surv_data = "../../metadata/COAD/COAD_cBioPortal_disease_free_survival_firehose.txt",
        survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
        categorical_covariates = c("history_colon_polyps"),
        values_not_considered = c("unknown"),
        timerange_cat = list(c(0, 5)),
        numeric_covariates = c("PC4", "age", "mutation_burden"),
        timerange_cont = list(c(0, 5), c(0, 5), c(0, 5))
    ), 
    output_file = "../../results/survival_analysis/COAD_selectedTumor_COX_DFS_PC4.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## Kaplan-Meyer

#---------------------------------------------------------------------------------------------------------------------------

## COAD noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/survival_analysis/km_analysis.Rmd", 
    params=list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt", 
                        "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_pca_tab.txt"),
        join = "columns",
        taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        surv_data = "../../metadata/COAD/COAD_cBioPortal_disease_free_survival_firehose.txt",
        survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
        categorical_covariates = c("history_colon_polyps"),
        values_not_considered = c("unknown"),
        timerange_cat = list(c(0, 5)),
        numeric_covariates = c("PC4", "age", "mutation_burden"),
        timerange_cont = list(c(0, 5), c(0, 5), c(0, 5)),
        break_val = 1
    ), 
    output_file = "../../results/survival_analysis/COAD_selectedTumor_KM_DFS_PC4.html"
)

rm(list=ls())
gc(full=TRUE)
