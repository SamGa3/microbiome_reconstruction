#---------------------------------------------------------------------------------------------------------------------------

## Taxa and clinical or molecular properties
## no batch correction

#---------------------------------------------------------------------------------------------------------------------------

# COAD noFFPE OnlyPrimary NoRiboZeroG 
# taxa prevalence 0.1

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/identification_related_species/taxa_compositions_wilc.Rmd", params=list(
  metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"),
  join_by="columns",
  taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
  cat_properties = c("gender", "bmi", "stage", "CMS", "history_of_other_malignancy", "side", "MSI_status", "CIMP_status", 
                    "history_colon_polyps"),
  values_not_considered = list("unknown", "unknown", "unknown", "unknown", c("unknown","inconsistency"), "unknown", 
                                "unknown", "unknown", "unknown"),
  cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score"),
  heatmap = TRUE,
  taxa_selection = "../../results/filters/Presence_more0.1samples_COAD.txt",
  taxa_selection_approach="intersect",
  table_path = "../../results/identification_related_species/tables/COAD_prevalence0.1_selectedTumor"
), output_file = "../../results/identification_related_species/COAD_prevalence0.1_selectedTumor_bacteria_species_compositions.html")

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# GBM OnlyPrimaryTumor NoDupl
# taxa prevalence 0.1

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/identification_related_species/taxa_compositions_wilc.Rmd", params=list(
  metadata = c("../../metadata/GBM/GBM_clinical_metadata.txt", "../../metadata/GBM/GBM_technical_metadata.txt"),
  join_by="columns",
  taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
  cat_properties = c("gender", "MSI_status"),
  values_not_considered = list("unknown", "unknown"),
  cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score"),
  heatmap = TRUE,
  taxa_selection = "../../results/filters/Presence_more0.1samples_GBM.txt",
  taxa_selection_approach="intersect",
  table_path = "../../results/identification_related_species/tables/GBM_prevalence0.1_selectedTumor"
), output_file = "../../results/identification_related_speciesGBM_prevalence0.1_selectedTumor_bacteria_species_compositions.html")

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# LUAD noFFPE OnlyPrimary NoRiboZeroG 
# taxa prevalence 0.1

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/identification_related_species/taxa_compositions_wilc.Rmd", params=list(
  metadata = c("../../metadata/LUAD/LUAD_clinical_metadata.txt", "../../metadata/LUAD/LUAD_technical_metadata.txt"),
  join_by="columns",
  taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
  cat_properties = c("gender", "stage", "history_of_other_malignancy", "MSI_status"),
  values_not_considered = list("unknown", "unknown", c("unknown","inconsistency"), "unknown"),
  cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score"),
  heatmap = TRUE,
  taxa_selection = "../../results/filters/Presence_more0.1samples_LUAD.txt",
  taxa_selection_approach="intersect",
  table_path = "../../results/identification_related_species/tables/LUAD_prevalence0.1_selectedTumor"
), output_file = "../../results/identification_related_species/LUAD_prevalence0.1_selectedTumor_bacteria_species_compositions.html")

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# LUSC OnlyPrimary NoDupl 
# taxa prevalence 0.1

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/identification_related_species/taxa_compositions_wilc.Rmd", params=list(
  metadata = c("../../metadata/LUSC/LUSC_clinical_metadata.txt", "../../metadata/LUSC/LUSC_technical_metadata.txt"),
  join_by="columns",
  taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
  cat_properties = c("gender", "stage", "MSI_status"),
  values_not_considered = list("unknown", "unknown", "unknown"),
  cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score"),
  heatmap = TRUE,
  taxa_selection = "../../results/filters/Presence_more0.1samples_LUSC.txt",
  taxa_selection_approach="intersect",
  table_path = "../../results/identification_related_species/tables/LUSC_prevalence0.1_selectedTumor"
), output_file = "../../results/identification_related_species/LUSC_prevalence0.1_selectedTumor_bacteria_species_compositions.html")

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# HNSC noFFPE OnlyPrimary
# taxa prevalence 0.1

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/identification_related_species/taxa_compositions_wilc.Rmd", params=list(
  metadata = c("../../metadata/HNSC/HNSC_clinical_metadata.txt", "../../metadata/HNSC/HNSC_technical_metadata.txt"),
  join_by="columns",
  taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
  cat_properties = c("gender", "stage", "history_of_other_malignancy", "MSI_status"),
  values_not_considered = list("unknown", "unknown", c("unknown","inconsistency"), "unknown"),
  cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score"),
  heatmap = TRUE,
  taxa_selection = "../../results/filters/Presence_more0.1samples_HNSC.txt",
  taxa_selection_approach="intersect",
  table_path = "../../results/identification_related_species/tables/HNSC_prevalence0.1_selectedTumor"
), output_file = "../../results/identification_related_species/HNSC_prevalence0.1_selectedTumor_bacteria_species_compositions.html")

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# OV OnlyPrimary NomirVana
# taxa prevalence 0.1

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/identification_related_species/taxa_compositions_wilc.Rmd", params=list(
  metadata = c("../../metadata/OV/OV_clinical_metadata.txt", "../../metadata/OV/OV_technical_metadata.txt"),
  join_by="columns",
  taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
  cat_properties = c("stage", "MSI_status"),
  values_not_considered = list("unknown", "unknown"),
  cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score"),
  heatmap = TRUE,
  taxa_selection = "../../results/filters/Presence_more0.1samples_OV.txt",
  taxa_selection_approach="intersect",
  table_path = "../../results/identification_related_speciesables/OV_prevalence0.1_selectedTumor"
), output_file = "../../results/identification_related_speciese/OV_prevalence0.1_selectedTumor_bacteria_species_compositions.html")

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# READ OnlyPrimary
# taxa prevalence 0.1

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/identification_related_species/taxa_compositions_wilc.Rmd", params=list(
  metadata = c("../../metadata/READ/READ_clinical_metadata.txt", "../../metadata/READ/READ_technical_metadata.txt"),
  join_by="columns",
  taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
  cat_properties = c("gender", "bmi", "stage", "CMS", "history_of_other_malignancy", "MSI_status", "CIMP_status", 
                    "history_colon_polyps"),
  values_not_considered = list("unknown", "unknown", "unknown", "unknown", c("unknown","inconsistency"), "unknown", 
                                "unknown", "unknown"),
  cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score"),
  heatmap = TRUE,
  taxa_selection = "../../results/filters/Presence_more0.1samples_READ.txt",
  taxa_selection_approach="intersect",
  table_path = "../../results/identification_related_species/tables/READ_prevalence0.1_OnlyPrimary"
), output_file = "../../results/identification_related_species/READ_prevalence0.1_selectedTumor_bacteria_species_compositions.html")

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# SKCM OnlyPrimary 
# taxa prevalence 0.1

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/identification_related_species/taxa_compositions_wilc.Rmd", params=list(
  metadata = c("../../metadata/SKCM/SKCM_clinical_metadata.txt", "../../metadata/SKCM/SKCM_technical_metadata.txt"),
  join_by="columns",
  taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
  cat_properties = c("gender", "stage", "history_of_other_malignancy"),
  values_not_considered = list("unknown", c("unknown", "NOS"), c("unknown","inconsistency")),
  cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score"),
  heatmap = TRUE,
  taxa_selection = "../../results/filters/Presence_more0.1samples_SKCM.txt",
  taxa_selection_approach="intersect",
  table_path = "../../results/identification_related_species/tables/SKCM_prevalence0.1_OnlyPrimary"
), output_file = "../../results/identification_related_species/SKCM_prevalence0.1_selectedTumor_bacteria_species_compositions.html")

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# BRCA noFFPE OnlyPrimary NoRiboZeroG 
# taxa prevalence 0.1

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/identification_related_species/taxa_compositions_wilc.Rmd", params=list(
  metadata = c("../../metadata/BRCA/BRCA_clinical_metadata.txt", "../../metadata/BRCA/BRCA_technical_metadata.txt"),
  join_by="columns",
  taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
  cat_properties = c("gender", "stage", "history_of_other_malignancy", "MSI_status"),
  values_not_considered = list("unknown", "unknown", c("unknown","inconsistency"), "unknown"),
  cont_properties = c("percent_normal_cells", "age_at_diagnosis", "mutation_burden", "stemness", "aneuploidy_score"),
  heatmap = TRUE,
  taxa_selection = "../../results/filters/Presence_more0.1samples_BRCA.txt",
  taxa_selection_approach="intersect",
  table_path = "../../results/identification_related_species/tables/BRCA_prevalence0.1_selectedTumor"
), output_file = "../../results/identification_related_species/BRCA_prevalence0.1_selectedTumor_bacteria_species_compositions.html")

rm(list=ls())
gc(full=TRUE)
