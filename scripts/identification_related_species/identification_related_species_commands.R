### COMMANDS

#---------------------------------------------------------------------------------------------------------------------------

## Taxa and clinical or molecular properties

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/identification_related_species/taxa_compositions.Rmd", 
    params=list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join_by = "columns",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        cat_properties = c("CMS", "side", "MSI_status"),
        cat_correction = rep("corr_plate_id", 3),
        cont_properties = c("Mast_cells_activated", "Mast_cells_resting"),
        cont_correction = rep("corr_plate_id", 2),
        values_not_considered = list(c("unknown", "NOLBL"), "unknown", "unknown"),
        taxa_selection = c("../../results/filters/Presence_more0.1samples_COAD.txt",
                            "../../results/filters/HighMeanVsRest_COAD_vs_GBM_LUAD_LUSC_HNSC_OV_SKCM.txt"),
        taxa_selection_approach="intersect",
        palette = c("CMS", "Left_right", "MSI", "jama", "nejm"),
        table_path = "../../results/identification_related_species/tables/COAD_selectedTumor_3filters"
    ), 
    output_file = "../../results/identification_related_species/COAD_selectTumor_3filters_bacteria_species_compositions.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## COX 

#---------------------------------------------------------------------------------------------------------------------------

# top 100 bacteria species used for the PC4

mat=read.csv("results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_rotations.txt", sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
bact=paste("X", mat$PC4, sep="")[1:100]

rmarkdown::render("scripts/survival_analysis/cox_univariate_analysis_on_microbes_tab.Rmd", 
  params=list(
    taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    surv_data = "../../metadata/COAD/COAD_cBioPortal_disease_free_survival_firehose.txt",
    survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
    numeric_covariates = bact,
    timerange_cont = rep(list(c(0, 5)), 100),
    output = "../../results/identification_related_species/tables/COAD_selectedTumor_top100PC4_bacteria_species_DFS_"
  ), 
  output_file = "../../results/identification_related_species/COAD_selectedTumor_cox_univa_top100PC4.html"
)

rm(list=ls())
gc(full=TRUE)
