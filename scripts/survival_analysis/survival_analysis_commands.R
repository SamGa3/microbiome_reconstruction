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
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/COAD_selectedTumor_COX_OS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/COAD_selectedTumor_COX_OS_PC123456.html"
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
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/COAD_selectedTumor_COX_DFS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/COAD_selectedTumor_COX_DFS_PC123456.html"
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

## GBM OnlyPrimaryTumor NoDupl

#---------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
gc(full=TRUE)

# PC1 2 3 4 5 6
rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/GBM/GBM_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/raw/merged_unamb_score/GBM/tables/GBM_selectedTumor_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../bacteria/RNAseq/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/GBM/GBM_cBioPortal_overall_survival_firehose.txt",
    survival_analysis = c("OS_YEARS", "patient_status", "OS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/GBM_OnlyPrimaryNoDupl_COX_OS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/GBM_selectedTumor_COX_OS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/GBM/GBM_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/raw/merged_unamb_score/GBM/tables/GBM_selectedTumor_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../bacteria/RNAseq/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/GBM/GBM_cBioPortal_disease_free_survival_firehose.txt",
    survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/GBM_OnlyPrimaryNoDupl_COX_DFS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/GBM_selectedTumor_COX_DFS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## LUAD noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

# PC1 2 3 4 5 6
rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/LUAD/LUAD_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/tables/LUAD_ComBat_corr_plate_id_selectedTumor_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/LUAD/LUAD_cBioPortal_overall_survival_firehose.txt",
    survival_analysis = c("OS_YEARS", "patient_status", "OS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/LUAD_selectedTumor_COX_OS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/LUAD_selectedTumor_COX_OS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/LUAD/LUAD_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/tables/LUAD_ComBat_corr_plate_id_selectedTumor_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/LUAD/LUAD_cBioPortal_disease_free_survival_firehose.txt",
    survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/LUAD_selectedTumor_COX_DFS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/LUAD_selectedTumor_COX_DFS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## LUSC OnlyPrimary NoDupl

#---------------------------------------------------------------------------------------------------------------------------

# PC1 2 3 4 5 6
rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/LUSC/LUSC_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/tables/LUSC_ComBat_corr_plate_id_OnlyPrimaryNoDupl_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUSC_ComBat_corr_plate_id_OnlyPrimaryNoDupl_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/LUSC/LUSC_cBioPortal_overall_survival_firehose.txt",
    survival_analysis = c("OS_YEARS", "patient_status", "OS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/LUSC_OnlyPrimaryNoDupl_COX_OS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/LUSC_selectedTumor_COX_OS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/LUSC/LUSC_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/tables/LUSC_ComBat_corr_plate_id_OnlyPrimaryNoDupl_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUSC_ComBat_corr_plate_id_OnlyPrimaryNoDupl_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/LUSC/LUSC_cBioPortal_disease_free_survival_firehose.txt",
    survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/LUSC_OnlyPrimaryNoDupl_COX_DFS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/LUSC_selectedTumor_COX_DFS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## HNSC noFFPE OnlyPrimary

#---------------------------------------------------------------------------------------------------------------------------

# PC1 2 3 4 5 6
rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/HNSC/HNSC_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/tables/HNSC_ComBat_corr_plate_id_noFFPEOnlyPrimary_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/HNSC_ComBat_corr_plate_id_noFFPEOnlyPrimary_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/HNSC/HNSC_cBioPortal_overall_survival_firehose.txt",
    survival_analysis = c("OS_YEARS", "patient_status", "OS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/HNSC_noFFPEOnlyPrimary_COX_OS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/HNSC_selectedTumor_COX_OS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/HNSC/HNSC_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/tables/HNSC_ComBat_corr_plate_id_noFFPEOnlyPrimary_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/HNSC_ComBat_corr_plate_id_noFFPEOnlyPrimary_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/HNSC/HNSC_cBioPortal_disease_free_survival_firehose.txt",
    survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/HNSC_noFFPEOnlyPrimary_COX_DFS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/HNSC_selectedTumor_COX_DFS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## OV noFFPE OnlyPrimary

#---------------------------------------------------------------------------------------------------------------------------

# PC1 2 3 4 5 6
rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/OV/OV_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/tables/OV_ComBat_corr_plate_id_OnlyPrimaryNomirVana_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/OV_ComBat_corr_plate_id_OnlyPrimaryNomirVana_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/OV/OV_cBioPortal_overall_survival_firehose.txt",
    survival_analysis = c("OS_YEARS", "patient_status", "OS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/OV_OnlyPrimaryNomirVana_COX_OS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/OV_selectedTumor_COX_OS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/OV/OV_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/tables/OV_ComBat_corr_plate_id_OnlyPrimaryNomirVana_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/OV_ComBat_corr_plate_id_OnlyPrimaryNomirVana_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/OV/OV_cBioPortal_disease_free_survival_firehose.txt",
    survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/OV_OnlyPrimaryNomirVana_COX_DFS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/OV_selectedTumor_COX_DFS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## SKCM OnlyPrimary

#---------------------------------------------------------------------------------------------------------------------------

# PC1 2 3 4 5 6
rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/SKCM/SKCM_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/tables/SKCM_ComBat_corr_plate_id_OnlyPrimary_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/SKCM_ComBat_corr_plate_id_OnlyPrimary_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/SKCM/SKCM_cBioPortal_overall_survival_firehose.txt",
    survival_analysis = c("OS_YEARS", "patient_status", "OS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/SKCM_OnlyPrimary_COX_OS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/SKCM_selectedTumor_COX_OS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/SKCM/SKCM_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/tables/SKCM_ComBat_corr_plate_id_OnlyPrimary_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/SKCM_ComBat_corr_plate_id_OnlyPrimary_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/SKCM/SKCM_cBioPortal_disease_free_survival_firehose.txt",
    survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/SKCM_OnlyPrimary_COX_DFS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/SKCM_selectedTumor_COX_DFS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## READ OnlyPrimary

#---------------------------------------------------------------------------------------------------------------------------

# PC1 2 3 4 5 6
rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/READ/READ_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/tables/READ_ComBat_corr_plate_id_OnlyPrimary_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_OnlyPrimary_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/READ/READ_cBioPortal_overall_survival_firehose.txt",
    survival_analysis = c("OS_YEARS", "patient_status", "OS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/READ_OnlyPrimary_COX_OS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/READ_selectedTumor_COX_OS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/READ/READ_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/tables/READ_ComBat_corr_plate_id_OnlyPrimary_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_OnlyPrimary_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/READ/READ_cBioPortal_disease_free_survival_firehose.txt",
    survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/READ_OnlyPrimary_COX_DFS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/READ_selectedTumor_COX_DFS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## BRCA noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

# PC1 2 3 4 5 6
rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/BRCA/BRCA_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/BRCA/tables/BRCA_ComBat_corr_plate_id_selectedTumor_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/BRCA_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/BRCA/BRCA_cBioPortal_overall_survival_firehose.txt",
    survival_analysis = c("OS_YEARS", "patient_status", "OS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/BRCA_selectedTumor_COX_OS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/BRCA_selectedTumor_COX_OS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
  params=list(
    metadata = c("../../metadata/BRCA/BRCA_clinical_metadata.txt",
                  "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/BRCA/tables/BRCA_ComBat_corr_plate_id_selectedTumor_pca_tab.txt"),
    join = c("columns"),
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/BRCA_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/BRCA/BRCA_cBioPortal_disease_free_survival_firehose.txt",
    survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
    categorical_covariates = "",
    values_not_considered = "",
    timerange_cat = "",
    numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
    timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5)),
    table_path = "../../results/survival_analysis/tables/BRCA_selectedTumor_COX_DFS_PC123456_"
  ), 
  output_file = "../../results/survival_analysis/BRCA_selectedTumor_COX_DFS_PC123456.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# Correct p values

rm(list=ls())
gc(full=TRUE)

COAD_OS=read.csv("../../results/survival_analysis/COAD/tables/COAD_selectedTumor_COX_OS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
COAD_OS=data.frame(cancer_type=rep("COAD", nrow(COAD_OS)), surv=rep("OS", nrow(COAD_OS)), PCs=rownames(COAD_OS), COAD_OS)
COAD_DFS=read.csv("../../results/survival_analysis/COAD/tables/COAD_selectedTumor_COX_DFS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
COAD_DFS=data.frame(cancer_type=rep("COAD", nrow(COAD_DFS)), surv=rep("DFS", nrow(COAD_OS)), PCs=rownames(COAD_DFS), COAD_DFS)

GBM_OS=read.csv("../../results/survival_analysis/GBM/tables/GBM_OnlyPrimaryNoDupl_COX_OS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
GBM_OS=data.frame(cancer_type=rep("GBM", nrow(GBM_OS)), surv=rep("OS", nrow(GBM_OS)), PCs=rownames(GBM_OS), GBM_OS)
GBM_DFS=read.csv("../../results/survival_analysis/GBM/tables/GBM_OnlyPrimaryNoDupl_COX_DFS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
GBM_DFS=data.frame(cancer_type=rep("GBM", nrow(GBM_DFS)), surv=rep("DFS", nrow(GBM_OS)), PCs=rownames(GBM_DFS), GBM_DFS)

HNSC_OS=read.csv("../../results/survival_analysis/HNSC/tables/HNSC_noFFPEOnlyPrimary_COX_OS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
HNSC_OS=data.frame(cancer_type=rep("HNSC", nrow(HNSC_OS)), surv=rep("OS", nrow(HNSC_OS)), PCs=rownames(HNSC_OS), HNSC_OS)
HNSC_DFS=read.csv("../../results/survival_analysis/HNSC/tables/HNSC_noFFPEOnlyPrimary_COX_DFS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
HNSC_DFS=data.frame(cancer_type=rep("HNSC", nrow(HNSC_DFS)), surv=rep("DFS", nrow(HNSC_OS)), PCs=rownames(HNSC_DFS), HNSC_DFS)

LUAD_OS=read.csv("../../results/survival_analysis/LUAD/tables/LUAD_selectedTumor_COX_OS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
LUAD_OS=data.frame(cancer_type=rep("LUAD", nrow(LUAD_OS)), surv=rep("OS", nrow(LUAD_OS)), PCs=rownames(LUAD_OS), LUAD_OS)
LUAD_DFS=read.csv("../../results/survival_analysis/LUAD/tables/LUAD_selectedTumor_COX_DFS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
LUAD_DFS=data.frame(cancer_type=rep("LUAD", nrow(LUAD_DFS)), surv=rep("DFS", nrow(LUAD_OS)), PCs=rownames(LUAD_DFS), LUAD_DFS)

LUSC_OS=read.csv("../../results/survival_analysis/LUSC/tables/LUSC_OnlyPrimaryNoDupl_COX_OS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
LUSC_OS=data.frame(cancer_type=rep("LUSC", nrow(LUSC_OS)), surv=rep("OS", nrow(LUSC_OS)), PCs=rownames(LUSC_OS), LUSC_OS)
LUSC_DFS=read.csv("../../results/survival_analysis/LUSC/tables/LUSC_OnlyPrimaryNoDupl_COX_DFS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
LUSC_DFS=data.frame(cancer_type=rep("LUSC", nrow(LUSC_DFS)), surv=rep("DFS", nrow(LUSC_OS)), PCs=rownames(LUSC_DFS), LUSC_DFS)

OV_OS=read.csv("../../results/survival_analysis/OV/tables/OV_OnlyPrimaryNomirVana_COX_OS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
OV_OS=data.frame(cancer_type=rep("OV", nrow(OV_OS)), surv=rep("OS", nrow(OV_OS)), PCs=rownames(OV_OS), OV_OS)
OV_DFS=read.csv("../../results/survival_analysis/OV/tables/OV_OnlyPrimaryNomirVana_COX_DFS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
OV_DFS=data.frame(cancer_type=rep("OV", nrow(OV_DFS)), surv=rep("DFS", nrow(OV_OS)), PCs=rownames(OV_DFS), OV_DFS)

READ_OS=read.csv("../../results/survival_analysis/READ/tables/READ_OnlyPrimary_COX_OS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
READ_OS=data.frame(cancer_type=rep("READ", nrow(READ_OS)), surv=rep("OS", nrow(READ_OS)), PCs=rownames(READ_OS), READ_OS)
READ_DFS=read.csv("../../results/survival_analysis/READ/tables/READ_OnlyPrimary_COX_DFS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
READ_DFS=data.frame(cancer_type=rep("READ", nrow(READ_DFS)), surv=rep("DFS", nrow(READ_OS)), PCs=rownames(READ_DFS), READ_DFS)

SKCM_OS=read.csv("../../results/survival_analysis/SKCM/tables/SKCM_OnlyPrimary_COX_OS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
SKCM_OS=data.frame(cancer_type=rep("SKCM", nrow(SKCM_OS)), surv=rep("OS", nrow(SKCM_OS)), PCs=rownames(SKCM_OS), SKCM_OS)
SKCM_DFS=read.csv("../../results/survival_analysis/SKCM/tables/SKCM_OnlyPrimary_COX_DFS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
SKCM_DFS=data.frame(cancer_type=rep("SKCM", nrow(SKCM_DFS)), surv=rep("DFS", nrow(SKCM_OS)), PCs=rownames(SKCM_DFS), SKCM_DFS)

BRCA_OS=read.csv("../../results/survival_analysis/BRCA/tables/BRCA_OnlyPrimary_COX_OS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
BRCA_OS=data.frame(cancer_type=rep("BRCA", nrow(BRCA_OS)), surv=rep("OS", nrow(BRCA_OS)), PCs=rownames(BRCA_OS), BRCA_OS)
BRCA_DFS=read.csv("../../results/survival_analysis/BRCA/tables/BRCA_OnlyPrimary_COX_DFS_PC123456_univ.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
BRCA_DFS=data.frame(cancer_type=rep("BRCA", nrow(BRCA_DFS)), surv=rep("DFS", nrow(BRCA_OS)), PCs=rownames(BRCA_DFS), BRCA_DFS)

os_tab=do.call(rbind, list(COAD_OS, GBM_OS, HNSC_OS, LUAD_OS, LUSC_OS, OV_OS, READ_OS, SKCM_OS, BRCA_OS))
os_tab=data.frame(os_tab, q.value=p.adjust(os_tab[,"p.value"]))
write.table(os_tab, file="../../results/survival_analysis/total_cancer_types_OS.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

dfs_tab=do.call(rbind, list(COAD_DFS, GBM_DFS, HNSC_DFS, LUAD_DFS, LUSC_DFS, OV_DFS, READ_DFS, SKCM_DFS, BRCA_DFS))
dfs_tab=data.frame(dfs_tab, q.value=p.adjust(dfs_tab[,"p.value"]))
write.table(dfs_tab, file="../../results/survival_analysis/total_cancer_types_DFS.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

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
