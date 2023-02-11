
rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## ALL TISSUES (COAD GBM LUAD LUSC HNSC OV READ SKCM)

# COAD noFFPE OnlyPrimary NoRiboZeroG
# GBM OnlyPrimary NoDupl
# LUAD noFFPE OnlyPrimary NoRiboZeroG
# LUSC OnlyPrimary NoDupl
# HNSC noFFPE OnlyPrimary
# OV OnlyPrimary NomirVana
# READ OnlyPrimary
# SKCM OnlyPrimary
# BRCA noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/GBM/GBM_clinical_metadata.txt", 
                      "../../metadata/LUAD/LUAD_clinical_metadata.txt", "../../metadata/LUSC/LUSC_clinical_metadata.txt", 
                      "../../metadata/HNSC/HNSC_clinical_metadata.txt", "../../metadata/OV/OV_clinical_metadata.txt", 
                      "../../metadata/READ/READ_clinical_metadata.txt", "../../metadata/SKCM/SKCM_clinical_metadata.txt",
                      "../../metadata/BRCA/BRCA_clinical_metadata.txt"
                    ),
        join = "rows",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
                ),
        cat_properties = c("project_id"),
        values_not_considered = list("unknown"),
        cont_properties = c(""),
        picture_3d=list(list(path="../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/all/images/COAD_GBM_LUAD_LUSC_HNSC_OV_READ_SKCM_BRCA", feat="project_id", angle=20, PCs=c(1,2,3))),
        palette = c("locuszoom2")
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/all/COAD_GBM_LUAD_LUSC_HNSC_OV_READ_SKCM_BRCA_selectedTumor_property_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## ALL TISSUES (IEO + COAD GBM LUAD LUSC HNSC OV READ SKCM BRCA)

# IEO list OnlySelected
# COAD noFFPE NormalPrimary NoRiboZeroG
# GBM NormalPrimary NoDupl
# LUAD noFFPE NormalPrimary NoRiboZeroG
# LUSC NormalPrimary NoDupl
# HNSC noFFPE NormalPrimary
# OV OnlyPrimary NomirVana
# READ NormalPrimary
# SKCM NOrmalPrimary
# BRCA noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/IEO/IEO_clinical_metadata.txt",
                      "../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/GBM/GBM_clinical_metadata.txt", 
                      "../../metadata/LUAD/LUAD_clinical_metadata.txt", "../../metadata/LUSC/LUSC_clinical_metadata.txt", 
                      "../../metadata/HNSC/HNSC_clinical_metadata.txt", "../../metadata/OV/OV_clinical_metadata.txt", 
                      "../../metadata/READ/READ_clinical_metadata.txt", "../../metadata/SKCM/SKCM_clinical_metadata.txt",
                      "../../metadata/BRCA/BRCA_clinical_metadata.txt"
                    ),
        join = "rows",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
                ),
        cat_properties = c("project_id"),
        values_not_considered = list("unknown"),
        cont_properties = c(""),
        picture_3d=list(list(path="../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/all/images/IEO_COAD_GBM_LUAD_LUSC_HNSC_OV_READ_SKCM_BRCA", feat="project_id", angle=20, PCs=c(1,2,3))),
        palette = c("COAD_IEO_rest")
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/all/IEO_COAD_GBM_LUAD_LUSC_HNSC_OV_READ_SKCM_BRCA_selectedTumor_property_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## COAD

#---------------------------------------------------------------------------------------------------------------------------

# OnlyNormal
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedNormal_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "bmi", "stage", "CMS", "history_of_other_malignancy", "side", "MSI_status", "CIMP_status", "history_colon_polyps"),
        values_not_considered = list("unknown", "unknown", "unknown", c("unknown", "NOLBL"), c("unknown","inconsistency"), "unknown", "unknown", "unknown", "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/COAD/COAD_selectedNormal_property_association.html"
)

rm(list=ls())
gc(full=TRUE)

# noFFPE OnlyPrimary NoRiboZeroG
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("corr_plate_id", "read_length"),
        values_not_considered = list("unknown", "unknown"),
        cont_properties = c(""),
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        pca2prop = list(c(cat="corr_plate_id"), c(cat="read_length"), c("brewer_Paired_Dark2_Set2")),
        palette = c("default", "default")
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_property_association.html"
)

rm(list=ls())
gc(full=TRUE)

# OnlyCoupled noFFPE OnlyPrimary NoRiboZeroG 
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedCoupledTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "bmi", "stage", "CMS", "history_of_other_malignancy", "side", "MSI_status", "CIMP_status", "history_colon_polyps"),
        values_not_considered = list("unknown", "unknown", "unknown", c("unknown", "NOLBL"), c("unknown","inconsistency"), "unknown", "unknown", "unknown", "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/COAD/COAD_selectedCoupledTumor_property_association.html"
)

rm(list=ls())
gc(full=TRUE)

# noFFPE OnlyPrimary NoRiboZeroG ReadLength48 
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumorReadLength48_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "bmi", "stage", "CMS", "history_of_other_malignancy", "side", "MSI_status", "CIMP_status", "history_colon_polyps"),
        values_not_considered = list("unknown", "unknown", "unknown", c("unknown", "NOLBL"), c("unknown","inconsistency"), "unknown", "unknown", "unknown", "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/COAD/COAD_selectedTumorReadLength48_property_association.html"
)

rm(list=ls())
gc(full=TRUE)

# noFFPE OnlyPrimary NoRiboZeroG ReadLength76
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumorReadLength76_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "bmi", "stage", "CMS", "history_of_other_malignancy", "side", "MSI_status", "CIMP_status", "history_colon_polyps"),
        values_not_considered = list("unknown", "unknown", "unknown", c("unknown", "NOLBL"), c("unknown","inconsistency"), "unknown", "unknown", "unknown", "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/COAD/COAD_selectedTumorReadLength76_property_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## GBM OnlyPrimaryTumor NoDupl

#---------------------------------------------------------------------------------------------------------------------------

# OnlyPrimaryTumor NoDupl
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/GBM/GBM_technical_metadata.txt", "../../metadata/GBM/GBM_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "MSI_status"),
        values_not_considered = list("unknown", "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_property_association.html"
)

rm(list=ls())
gc(full=TRUE)

# cibersort relative
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/GBM/GBM_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = "",
        values_not_considered = "",
        cont_properties = c("B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", 
                            "T_cells_follicular_helper", "T_cells_regulatory_Tregs", "T_cells_gamma_delta", "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", 
                            "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated", 
                            "Eosinophils", "Neutrophils")
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_immuneInfiltrationRelative_pbelow05_association.html"
)

rm(list=ls())
gc(full=TRUE)

# mutations
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/GBM/GBM_mutation_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("PTEN", "TP53", "EGFR", "PIK3R1", "PIK3CA", "NF1", "RB1", "IDH1", "STAG2", "RPL5", "SLC26A3", "MUC17", "KEL", "CHD8", "NUP210L", 
                            "BRAF", "MAP3K1", "QKI", "AZGP1", "SETD2", #"CD1D", 
                            "DDX5", "IL4R", "BCOR", "PAN3", "AFM", "CD209", "SMC1A", "PTPN11", "GLT8D2", "EZR", "ODF4", "ZDHHC4", "PRKCD"
                          ),
        values_not_considered = rep("unknown", 33),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/GBM/tables/GBM_selectedTumor_mutation_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_mutation_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy whole chr
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/GBM/GBM_wholeChrAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
                          ),
        values_not_considered = rep("unknown", 22),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/GBM/tables/GBM_selectedTumor_wholeChrAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_wholeChrAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy armlevel
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/GBM/GBM_armAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                            "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", 
                            "chr13q", "chr14q", "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", 
                            "chr20p", "chr20q", "chr21q", "chr22q"
                          ),
        values_not_considered = rep("unknown", 39),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/GBM/tables/GBM_selectedTumor_armAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_armAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## COAD ComBat read_length

#---------------------------------------------------------------------------------------------------------------------------

# noFFPE OnlyPrimary NoRiboZeroG
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_read_length/merged_unamb_score_norm/COAD_ComBat_read_length_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "bmi", "stage", "CMS", "history_of_other_malignancy", "side", "MSI_status", "CIMP_status", "history_colon_polyps"),
        values_not_considered = list("unknown", "unknown", "unknown", c("unknown", "NOLBL"), c("unknown","inconsistency"), "unknown", "unknown", "unknown", "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score"),
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        pca2prop = list(c(cat="corr_plate_id"), c(cat="read_length"), c("brewer_Paired_Dark2_Set2"))
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_read_length/merged_unamb_score_norm/COAD/COAD_ComBatReadLength_selectedTumor_property_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## COAD ComBat corr_plate_id

#---------------------------------------------------------------------------------------------------------------------------

# noFFPE OnlyPrimary NoRiboZeroG
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "bmi", "stage", "CMS", "history_of_other_malignancy", "side", "MSI_status", "CIMP_status", "history_colon_polyps"),
        values_not_considered = list("unknown", "unknown", "unknown", c("unknown", "NOLBL"), c("unknown","inconsistency"), "unknown", "unknown", "unknown", "unknown"),
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score"),
        rotations_path="../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_rotations.txt",
        pca_matrix_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_pca_tab.txt",
        pca2prop = list(c(cat="corr_plate_id"), c(cat="read_length"), c("brewer_Paired_Dark2_Set2")),
        palette = c("default", "default", "default", "CMS", rep("default", 11))
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/COAD_selectedTumor_property_association.html"
)

rm(list=ls())
gc(full=TRUE)

# cibersort relative
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = "",
        values_not_considered = "",
        cont_properties = c("B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", 
                            "T_cells_follicular_helper", "T_cells_regulatory_Tregs", "T_cells_gamma_delta", "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", 
                            "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated", 
                            "Eosinophils", "Neutrophils")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/COAD_selectedTumor_immuneInfiltrationRelative_pbelow05_association.html"
)

rm(list=ls())
gc(full=TRUE)

# cibersort absolute
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_immuneInfiltrationAbsolute_pbelow05_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = "",
        values_not_considered = "",
        cont_properties = c("B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", 
                            "T_cells_follicular_helper", "T_cells_regulatory_Tregs", "T_cells_gamma_delta", "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", 
                            "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated", 
                            "Eosinophils", "Neutrophils")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/COAD_selectedTumor_immuneInfiltrationAbsolute_pbelow05_association.html"
)

rm(list=ls())
gc(full=TRUE)

# mutations
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_mutation_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("APC", "TP53", "KRAS", "PIK3CA", "FBXW7", "SMAD4", "TCF7L2", "NRAS", "BRAF", "SMAD2", "PCBP1",
                            "ARID1A", "ACVR1B", "ERBB3", "AXIN2", "CDC27", "CASP8", "ELF3", "NTN4", "B2M", "GOT1", "TRIM23", "SIRT4",
                            #"FAM123B", "KIAA1804", 
                            "ARID2", "MSH6", "MIER3", "MAP2K4", "CTNNB1", "ERBB2", "BCLAF1", "BCOR", "NR4A2", 
                            "RBM10", "IDH2", "PTEN", "CNBD1", "TRAF3", "CD70"
                          ),
        values_not_considered = rep("unknown", 38),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_selectedTumor_mutation_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/COAD_selectedTumor_mutation_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy whole chr
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_wholeChrAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
                          ),
        values_not_considered = rep("unknown", 22),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_selectedTumor_wholeChrAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/COAD_selectedTumor_wholeChrAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy armlevel
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_armAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                            "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", 
                            "chr13q", "chr14q", "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", 
                            "chr20p", "chr20q", "chr21q", "chr22q"
                          ),
        values_not_considered = rep("unknown", 39),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_selectedTumor_armAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/COAD_selectedTumor_armAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## LUAD ComBat corr_plate_id

#---------------------------------------------------------------------------------------------------------------------------

# noFFPE OnlyPrimary NoRiboZeroG
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/LUAD/LUAD_technical_metadata.txt", "../../metadata/LUAD/LUAD_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "stage", "history_of_other_malignancy", "MSI_status"),
        values_not_considered = list("unknown", "unknown", c("unknown","inconsistency"), "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_property_association.html"
)

# cibersort relative
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/LUAD/LUAD_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = "",
        values_not_considered = "",
        cont_properties = c("B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", 
                            "T_cells_follicular_helper", "T_cells_regulatory_Tregs", "T_cells_gamma_delta", "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", 
                            "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated", 
                            "Eosinophils", "Neutrophils")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_immuneInfiltrationRelative_pbelow05_association.html"
)

rm(list=ls())
gc(full=TRUE)

# mutations
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/LUAD/LUAD_mutation_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("TP53", "KRAS", "KEAP1", "EGFR", "STK11", "NF1", "SMARCA4", "BRAF", "CDKN2A", "MET", "RBM10", "PIK3CA", "RIT1", "U2AF1", 
                            "ATM", "ARID1A", "SLC4A5", "RB1", "ERBB2", "MAP2K1", "STX2", "NBPF1", #"MLL3", 
                            "FAT1", "APC", "ARID2", "SLC39A6", 
                            "ARHGAP35", "SMAD4", "CTNNB1", "CDK12", "MBD1", "FBXO16", "PBLD", "ZNF774", "NRAS", "CDKN1B"
                          ),
        values_not_considered = rep("unknown", 36),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/tables/LUAD_selectedTumor_mutation_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_mutation_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy whole chr
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/LUAD/LUAD_wholeChrAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
                          ),
        values_not_considered = rep("unknown", 22),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/tables/LUAD_selectedTumor_wholeChrAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_wholeChrAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy armlevel
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/LUAD/LUAD_armAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                            "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", 
                            "chr13q", "chr14q", "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", 
                            "chr20p", "chr20q", "chr21q", "chr22q"
                          ),
        values_not_considered = rep("unknown", 39),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/tables/LUAD_selectedTumor_armAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_armAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## LUSC ComBat corr_plate_id

#---------------------------------------------------------------------------------------------------------------------------

# OnlyPrimary NoDupl
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/LUSC/LUSC_technical_metadata.txt", "../../metadata/LUSC/LUSC_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "stage", "MSI_status"),
        values_not_considered = list("unknown", "unknown", "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_property_association.html"
)

# cibersort relative
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/LUSC/LUSC_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = "",
        values_not_considered = "",
        cont_properties = c("B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", 
                            "T_cells_follicular_helper", "T_cells_regulatory_Tregs", "T_cells_gamma_delta", "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", 
                            "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated", 
                            "Eosinophils", "Neutrophils")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_immuneInfiltrationRelative_pbelow05_association.html"
)

rm(list=ls())
gc(full=TRUE)

# mutations
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/LUSC/LUSC_mutation_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("TP53", #"MLL2", 
                            "CDKN2A", "PIK3CA", "NFE2L2", "KEAP1", "ARID1A", "RB1", "FBXW7", "HLA.A", "HRAS", "FAT1", 
                            "NF1", "PTEN", "NOTCH1", "ARHGAP35", "NSD1", "ASXL1", "EP300", "RASA1", "EGFR", "FUBP1", "FGFR3", "STK11", 
                            "SERPINB13"
                          ),
        values_not_considered = rep("unknown", 24),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/tables/LUSC_selectedTumor_mutation_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_mutation_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy whole chr
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/LUSC/LUSC_wholeChrAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
                          ),
        values_not_considered = rep("unknown", 22),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/tables/LUSC_selectedTumor_wholeChrAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_wholeChrAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy armlevel
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/LUSC/LUSC_armAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                            "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", 
                            "chr13q", "chr14q", "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", 
                            "chr20p", "chr20q", "chr21q", "chr22q"
                          ),
        values_not_considered = rep("unknown", 39),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/tables/LUSC_selectedTumor_armAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_armAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## HNSC ComBat corr_plate_id

#---------------------------------------------------------------------------------------------------------------------------

# noFFPE OnlyPrimary

# OnlyPrimary NoDupl
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/HNSC/HNSC_technical_metadata.txt", "../../metadata/HNSC/HNSC_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/HNSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "stage", "history_of_other_malignancy", "MSI_status"),
        values_not_considered = list("unknown", "unknown", c("unknown","inconsistency"), "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_property_association.html"
)

# cibersort relative
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/HNSC/HNSC_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/HNSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = "",
        values_not_considered = "",
        cont_properties = c("B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", 
                            "T_cells_follicular_helper", "T_cells_regulatory_Tregs", "T_cells_gamma_delta", "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", 
                            "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated", 
                            "Eosinophils", "Neutrophils")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_immuneInfiltrationRelative_pbelow05_association.html"
)

rm(list=ls())
gc(full=TRUE)

# mutations
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/HNSC/HNSC_mutation_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/HNSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("TP53", "CDKN2A", "FAT1", "PIK3CA", "NOTCH1", #"MLL2", 
                            "NSD1", "CASP8", "NFE2L2", "AJUBA", #"EPHA", 
                            "HRAS", "ZNF750",
                            "TGFBR2", "HLA.A", "PTEN", "RAC1", "RHOA", "EP300", "RASA1", "CTCF", "HLA.B", "B2M", "MAP4K3", "IPO7", "OTUD7A",
                            #"MLL3", 
                            "NCOR1", "ARHGAP35", "LCP1", "RB1", "ARID2", "ASXL1", "PBRM1", "KDM6A", "IRF6", "SMAD4", "RBM5", "ELF4",
                            "MAPK1"
                          ),
        values_not_considered = rep("unknown", 37),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/tables/HNSC_selectedTumor_mutation_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_mutation_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy whole chr
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/HNSC/HNSC_wholeChrAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/HNSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
                          ),
        values_not_considered = rep("unknown", 22),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/tables/HNSC_selectedTumor_wholeChrAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_wholeChrAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy armlevel
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/HNSC/HNSC_armAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/HNSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                            "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", 
                            "chr13q", "chr14q", "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", 
                            "chr20p", "chr20q", "chr21q", "chr22q"
                          ),
        values_not_considered = rep("unknown", 39),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/tables/HNSC_selectedTumor_armAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_armAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## OV ComBat corr_plate_id

#---------------------------------------------------------------------------------------------------------------------------

# OnlyPrimary NomirVana 
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/OV/OV_technical_metadata.txt", "../../metadata/OV/OV_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/OV_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("stage", "MSI_status"),
        values_not_considered = list("unknown", "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/OV_selectedTumor_property_association.html"
)

# cibersort relative
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/OV/OV_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/OV_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = "",
        values_not_considered = "",
        cont_properties = c("B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", 
                            "T_cells_follicular_helper", "T_cells_regulatory_Tregs", "T_cells_gamma_delta", "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", 
                            "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated", 
                            "Eosinophils", "Neutrophils")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/OV_selectedTumor_immuneInfiltrationRelative_pbelow05_association.html"
)

rm(list=ls())
gc(full=TRUE)

# mutations
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/OV/OV_mutation_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/OV_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("TP53", "BRCA1", "RB1", "CDK12", 
                            "NF1", 
                            "BRCA2", "COL5A3", "CREBBP", "SMARCB1", "ITGB7", "ERBB2", "ACVR2B", #"PDAP1", 
                            "KRAS"#, "NRAS"
                          ),
        values_not_considered = rep("unknown", 13),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/tables/OV_selectedTumor_mutation_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/OV_selectedTumor_mutation_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy whole chr
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/OV/OV_wholeChrAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/OV_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
                          ),
        values_not_considered = rep("unknown", 22),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/tables/OV_selectedTumor_wholeChrAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/OV_selectedTumor_wholeChrAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy armlevel
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/OV/OV_armAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/OV_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                            "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", 
                            "chr13q", "chr14q", "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", 
                            "chr20p", "chr20q", "chr21q", "chr22q"
                          ),
        values_not_considered = rep("unknown", 39),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/tables/OV_selectedTumor_armAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/OV/OV_selectedTumor_armAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## READ ComBat corr_plate_id

#---------------------------------------------------------------------------------------------------------------------------

# OnlyPrimary
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/READ/READ_technical_metadata.txt", "../../metadata/READ/READ_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "bmi", "stage", "CMS", "history_of_other_malignancy", "MSI_status", "CIMP_status", 
                              "history_colon_polyps"),
        values_not_considered = list("unknown", "unknown", "unknown", "unknown", c("unknown","inconsistency"), "unknown", 
                                       "unknown", "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/READ_selectedTumor_property_association.html"
)

# cibersort relative
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/READ/READ_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = "",
        values_not_considered = "",
        cont_properties = c("B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", 
                            "T_cells_follicular_helper", "T_cells_regulatory_Tregs", "T_cells_gamma_delta", "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", 
                            "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated", 
                            "Eosinophils", "Neutrophils")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/READ_selectedTumor_immuneInfiltrationRelative_pbelow05_association.html"
)

rm(list=ls())
gc(full=TRUE)

# mutations
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/READ/READ_mutation_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("APC", "TP53", "KRAS", "PIK3CA", "FBXW7", "SMAD4", "TCF7L2", "NRAS", "BRAF", "SMAD2", "PCBP1",
                            "ARID1A", "ACVR1B", "ERBB3", "AXIN2", "CDC27", "CASP8", "ELF3", "NTN4", "B2M", "GOT1", "TRIM23", "SIRT4",
                            #"FAM123B", "KIAA1804", 
                            "ARID2", "MSH6", "MIER3", "MAP2K4", "CTNNB1", "ERBB2", "BCLAF1", "BCOR", "NR4A2", 
                            "RBM10", "IDH2", "PTEN", "CNBD1", "TRAF3", "CD70"
                          ),
        values_not_considered = rep("unknown", 38),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/tables/READ_selectedTumor_mutation_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/READ_selectedTumor_mutation_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy whole chr
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/READ/READ_wholeChrAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
                          ),
        values_not_considered = rep("unknown", 22),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/tables/READ_selectedTumor_wholeChrAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/READ_selectedTumor_wholeChrAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy armlevel
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/READ/READ_armAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                            "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", 
                            "chr13q", "chr14q", "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", 
                            "chr20p", "chr20q", "chr21q", "chr22q"
                          ),
        values_not_considered = rep("unknown", 39),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/tables/READ_selectedTumor_armAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/READ/READ_selectedTumor_armAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## SKCM ComBat corr_plate_id

#---------------------------------------------------------------------------------------------------------------------------

# OnlyPrimary
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/SKCM/SKCM_technical_metadata.txt", "../../metadata/SKCM/SKCM_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/SKCM_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "stage", "history_of_other_malignancy"),
        values_not_considered = list("unknown", c("unknown", "NOS"), c("unknown","inconsistency")),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_property_association.html"
)

# cibersort relative
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/SKCM/SKCM_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/SKCM_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = "",
        values_not_considered = "",
        cont_properties = c("B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", 
                            "T_cells_follicular_helper", "T_cells_regulatory_Tregs", "T_cells_gamma_delta", "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", 
                            "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated", 
                            "Eosinophils", "Neutrophils")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_immuneInfiltrationRelative_pbelow05_association.html"
)

rm(list=ls())
gc(full=TRUE)

# mutations
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/SKCM/SKCM_mutation_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/SKCM_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("BRAF", "NRAS", "BCLAF1", "TP53", "CDKN2A", "RAC1",
                            "XIRP2", "MYOCD", "ALPK2", "PTEN", "PPP6C", "FRMD7", 
                            "OR4A16", "OR52N1", "WASF3", "CDK4", "LCTL", "STK19", "ACO1",
                            "ANK3", "MXRA5", "GPR116", "SLC12A8", "EGFL6", "PCDP1", 
                            "ARID2", #"CCDC148", 
                            "SNX31", "SLC44A4", "CTNNB1", #"MGC42105", 
                            "ZNF217", "RXRA", "NCOR1", "ALDH1B1", "ZNF750", "PIK3CA", "FDXACB1", "KIT"
                          ),
        values_not_considered = rep("unknown", 37),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/tables/SKCM_selectedTumor_mutation_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_mutation_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy whole chr
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/SKCM/SKCM_wholeChrAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/SKCM_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
                          ),
        values_not_considered = rep("unknown", 22),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/tables/SKCM_selectedTumor_wholeChrAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_wholeChrAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy armlevel
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/SKCM/SKCM_armAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/SKCM_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                            "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", 
                            "chr13q", "chr14q", "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", 
                            "chr20p", "chr20q", "chr21q", "chr22q"
                          ),
        values_not_considered = rep("unknown", 39),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/tables/SKCM_selectedTumor_armAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_armAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## BRCA ComBat corr_plate_id

#---------------------------------------------------------------------------------------------------------------------------

# OnlyPrimary
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/BRCA/BRCA_technical_metadata.txt", "../../metadata/BRCA/BRCA_clinical_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/BRCA_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("gender", "stage", "history_of_other_malignancy", "MSI_status"),
        values_not_considered = list("unknown", c("unknown", "NOS"), c("unknown","inconsistency")),
        cont_properties = c("percent_normal_cells", "age_at_diagnosis", "mutation_burden", "stemness", "aneuploidy_score")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_property_association.html"
)

# cibersort relative
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/BRCA/BRCA_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/BRCA_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = "",
        values_not_considered = "",
        cont_properties = c("B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", 
                            "T_cells_follicular_helper", "T_cells_regulatory_Tregs", "T_cells_gamma_delta", "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", 
                            "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated", 
                            "Eosinophils", "Neutrophils")
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_immuneInfiltrationRelative_pbelow05_association.html"
)

rm(list=ls())
gc(full=TRUE)

# mutations
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/BRCA/BRCA_mutation_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/BRCA_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("PIK3CA", "TP53", "GATA3", "MAP3K1", #"MLL3", 
                            "CDH1", "NCOR1", 
                            "MAP2K4", "PTEN", "RUNX1", "PIK3R1", "CTCF", "AKT1", "CBFB", 
                            "SPEN", "SF3B1", "ARID1A", "RB1", #"MLL", 
                            "KRAS", 
                            "TBX3", "ERBB2", "FOXA1", "MED23", "STAG2", "MYB", "TBL1XR1", "HIST1H3B", 
                            "CASP8", "CDKN1B", "CUL4B", "RAB40A", 
                            "ERBB3", "CDC42BPA", "SETDB1", "FGFR2", "GNPTAB", "EP300", "ACVR1B"
                            ),
        values_not_considered = rep("unknown", 37),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/BRCA/tables/BRCA_selectedTumor_mutation_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_mutation_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy whole chr
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/BRCA/BRCA_wholeChrAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/BRCA_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
                          ),
        values_not_considered = rep("unknown", 22),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/BRCA/tables/BRCA_selectedTumor_wholeChrAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_wholeChrAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

# aneuploidy armlevel
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/BRCA/BRCA_armAneuploidy_metadata.txt"),
        join = "columns",
        taxa = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/BRCA_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        cat_properties = c("chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q", "chr5p", "chr5q", "chr6p", "chr6q", 
                            "chr7p", "chr7q", "chr8p", "chr8q", "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q", 
                            "chr13q", "chr14q", "chr15q", "chr16p", "chr16q", "chr17p", "chr17q", "chr18p", "chr18q", "chr19p", "chr19q", 
                            "chr20p", "chr20q", "chr21q", "chr22q"
                          ),
        values_not_considered = rep("unknown", 39),
        cont_properties = c(""),
        heatmap_qtab_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/BRCA/tables/BRCA_selectedTumor_armAneuploidy_association.txt",
        heatmap_nPCs = 6
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_armAneuploidy_association.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## IEO ComBat corr_plate_id

#---------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
gc(full=TRUE)

# OnlySelected OnlyTumor
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/IEO/IEO_technical_metadata.txt", "../../metadata/IEO/IEO_clinical_metadata.txt"),
        taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        join = "columns",
        cat_properties = c("side"),
        values_not_considered = "unknown",
        cont_properties = "",
        palette = "Left_right"
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/IEO/IEO_selectedTumor_clinical_association.html"
)

rm(list=ls())
gc(full=TRUE)
