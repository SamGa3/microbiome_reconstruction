
rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## COAD noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/filters/filters.Rmd", 
    params = list(
        metadata_tissue = c(COAD="../../metadata/COAD/COAD_technical_metadata.txt"),
        metadata_comp = c(GBM="../../metadata/GBM/GBM_technical_metadata.txt", 
                            LUAD="../../metadata/LUAD/LUAD_technical_metadata.txt", 
                            LUSC="../../metadata/LUSC/LUSC_technical_metadata.txt", 
                            HNSC="../../metadata/HNSC/HNSC_technical_metadata.txt", 
                            OV="../../metadata/OV/OV_technical_metadata.txt",
                            SKCM="../../metadata/SKCM/SKCM_technical_metadata.txt",
                            BRCA="../../metadata/BRCA/BRCA_technical_metadata.txt"),
        taxa_tissue = c(COAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_comp = c(GBM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        HNSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        OV="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        SKCM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        BRCA="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        feat_tissue = "Project_ID",
        batch_feat = "plate_id",
        output_folder = "../../results/filters/",
        thr_presence = 0.1
    ), 
    output_file = "../../results/filters/COAD_selectedTumor_vs_GBM_LUAD_LUSC_HNSC_OV_SKCM_BRCA_selectedTumor_filters.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## GBM OnlyPrimary NoDupl

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/filters/filters.Rmd", 
    params = list(
        metadata_tissue = c(GBM="../../metadata/GBM/GBM_technical_metadata.txt"),
        metadata_comp = c(COAD="../../metadata/COAD/COAD_technical_metadata.txt", 
                            LUAD="../../metadata/LUAD/LUAD_technical_metadata.txt", 
                            LUSC="../../metadata/LUSC/LUSC_technical_metadata.txt", 
                            HNSC="../../metadata/HNSC/HNSC_technical_metadata.txt", 
                            OV="../../metadata/OV/OV_technical_metadata.txt",
                            READ="../../metadata/READ/READ_technical_metadata.txt",
                            SKCM="../../metadata/SKCM/SKCM_technical_metadata.txt",
                            BRCA="../../metadata/BRCA/BRCA_technical_metadata.txt"),
        taxa_tissue = c(GBM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_comp = c(COAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        HNSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        OV="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        SKCM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        READ="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        BRCA="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        feat_tissue = "Project_ID",
        batch_feat = "plate_id",
        output_folder = "../../results/filters/",
        thr_presence = 0.1
    ), 
    output_file = "../../results/filters/GBM_selectedTumor_vs_COAD_LUAD_LUSC_HNSC_OV_READ_SKCM_BRCA_selectedTumor_filters.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## LUAD noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/filters/filters.Rmd", 
    params = list(
        metadata_tissue = c(LUAD="../../metadata/LUAD/LUAD_technical_metadata.txt"),
        metadata_comp = c(COAD="../../metadata/COAD/COAD_technical_metadata.txt", 
                            GBM="../../metadata/GBM/GBM_technical_metadata.txt", 
                            HNSC="../../metadata/HNSC/HNSC_technical_metadata.txt", 
                            OV="../../metadata/OV/OV_technical_metadata.txt",
                            READ="../../metadata/READ/READ_technical_metadata.txt",
                            SKCM="../../metadata/SKCM/SKCM_technical_metadata.txt",
                            BRCA="../../metadata/BRCA/BRCA_technical_metadata.txt"),
        taxa_tissue = c(LUAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_comp = c(COAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        GBM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        HNSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        OV="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        READ="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        SKCM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        BRCA="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        feat_tissue = "Project_ID",
        batch_feat = "plate_id",
        output_folder = "../../results/filters/",
        thr_presence = 0.1
    ), 
    output_file = "../../results/filters/LUAD_selectedTumor_vs_COAD_GBM_HNSC_OV_READ_SKCM_BRCA_selectedTumor_filters.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## LUSC OnlyPrimary NoDupl

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/filters/filters.Rmd", 
    params = list(
        metadata_tissue = c(LUSC="../../metadata/LUSC/LUSC_technical_metadata.txt"),
        metadata_comp = c(COAD="../../metadata/COAD/COAD_technical_metadata.txt", 
                            GBM="../../metadata/GBM/GBM_technical_metadata.txt", 
                            HNSC="../../metadata/HNSC/HNSC_technical_metadata.txt", 
                            OV="../../metadata/OV/OV_technical_metadata.txt",
                            READ="../../metadata/READ/READ_technical_metadata.txt",
                            SKCM="../../metadata/SKCM/SKCM_technical_metadata.txt",
                            BRCA="../../metadata/BRCA/BRCA_technical_metadata.txt"),
        taxa_tissue = c(LUSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_comp = c(COAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        GBM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        HNSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        OV="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        READ="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        SKCM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        BRCA="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        feat_tissue = "Project_ID",
        batch_feat = "plate_id",
        output_folder = "../../results/filters/",
        thr_presence = 0.1
    ), 
    output_file = "../../results/filters/LUSC_selectedTumor_vs_COAD_GBM_HNSC_OV_SKCM_READ_BRCA_selectedTumor_filters.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## HNSC noFFPE OnlyPrimary

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/filters/filters.Rmd", 
    params = list(
        metadata_tissue = c(HNSC="../../metadata/HNSC/HNSC_technical_metadata.txt"),
        metadata_comp = c(COAD="../../metadata/COAD/COAD_technical_metadata.txt", 
                            GBM="../../metadata/GBM/GBM_technical_metadata.txt", 
                            LUAD="../../metadata/LUAD/LUAD_technical_metadata.txt", 
                            LUSC="../../metadata/LUSC/LUSC_technical_metadata.txt", 
                            OV="../../metadata/OV/OV_technical_metadata.txt",
                            READ="../../metadata/READ/READ_technical_metadata.txt",
                            SKCM="../../metadata/SKCM/SKCM_technical_metadata.txt",
                            BRCA="../../metadata/BRCA/BRCA_technical_metadata.txt"),
        taxa_tissue = c(HNSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_comp = c(COAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        GBM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        OV="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        READ="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        SKCM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        BRCA="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        feat_tissue = "Project_ID",
        batch_feat = "plate_id",
        output_folder = "../../results/filters/",
        thr_presence = 0.1
    ), 
    output_file = "../../results/filters/HNSC_selectedTumor_vs_COAD_GBM_LUAD_LUSC_OV_READ_SKCM_BRCA_selectedTumor_filters.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## OV OnlyPrimary NomirVana

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/filters/filters.Rmd", 
    params = list(
        metadata_tissue = c(OV="../../metadata/OV/OV_technical_metadata.txt"),
        metadata_comp = c(COAD="../../metadata/COAD/COAD_technical_metadata.txt", 
                            GBM="../../metadata/GBM/GBM_technical_metadata.txt", 
                            LUAD="../../metadata/LUAD/LUAD_technical_metadata.txt", 
                            LUSC="../../metadata/LUSC/LUSC_technical_metadata.txt", 
                            HNSC="../../metadata/HNSC/HNSC_technical_metadata.txt",
                            READ="../../metadata/READ/READ_technical_metadata.txt",
                            SKCM="../../metadata/SKCM/SKCM_technical_metadata.txt",
                            BRCA="../../metadata/BRCA/BRCA_technical_metadata.txt"),
        taxa_tissue = c(OV="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_comp = c(COAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        GBM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        HNSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        READ="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        SKCM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        BRCA="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        feat_tissue = "Project_ID",
        batch_feat = "plate_id",
        output_folder = "../../results/filters/",
        thr_presence = 0.1
    ), 
    output_file = "../../results/filters/OV_selectedTumor_vs_COAD_GBM_LUAD_LUSC_HNSC_READ_SKCM_BRCA_selectedTumor_filters.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## READ OnlyPrimary

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/filters/filters.Rmd", 
    params = list(
        metadata_tissue = c(READ="../../metadata/READ/READ_technical_metadata.txt"),
        metadata_comp = c(GBM="../../metadata/GBM/GBM_technical_metadata.txt", 
                            LUAD="../../metadata/LUAD/LUAD_technical_metadata.txt", 
                            LUSC="../../metadata/LUSC/LUSC_technical_metadata.txt", 
                            HNSC="../../metadata/HNSC/HNSC_technical_metadata.txt",
                            OV="../../metadata/OV/OV_technical_metadata.txt",
                            SKCM="../../metadata/SKCM/SKCM_technical_metadata.txt",
                            BRCA="../../metadata/BRCA/BRCA_technical_metadata.txt"),
        taxa_tissue = c(READ="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_comp = c(GBM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        HNSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        OV="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        SKCM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        BRCA="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        feat_tissue = "Project_ID",
        batch_feat = "plate_id",
        output_folder = "../../results/filters/",
        thr_presence = 0.1
    ), 
    output_file = "../../results/filters/READ_selectedTumor_vs_GBM_LUAD_LUSC_HNSC_OV_SKCM_BRCA_selectedTumor_filters.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## SKCM OnlyPrimary

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/filters/filters.Rmd", 
    params = list(
        metadata_tissue = c(SKCM="../../metadata/SKCM/SKCM_technical_metadata.txt"),
        metadata_comp = c(COAD="../../metadata/COAD/COAD_technical_metadata.txt", 
                            GBM="../../metadata/GBM/GBM_technical_metadata.txt", 
                            LUAD="../../metadata/LUAD/LUAD_technical_metadata.txt", 
                            LUSC="../../metadata/LUSC/LUSC_technical_metadata.txt", 
                            HNSC="../../metadata/HNSC/HNSC_technical_metadata.txt",
                            OV="../../metadata/OV/OV_technical_metadata.txt",
                            READ="../../metadata/READ/READ_technical_metadata.txt",
                            BRCA="../../metadata/BRCA/BRCA_technical_metadata.txt"),
        taxa_tissue = c(SKCM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_comp = c(COAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        GBM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        HNSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        OV="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        READ="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        BRCA="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        feat_tissue = "Project_ID",
        batch_feat = "plate_id",
        output_folder = "../../results/filters/",
        thr_presence = 0.1
    ), 
    output_file = "../../results/filters/SKCM_selectedTumor_vs_COAD_GBM_LUAD_LUSC_HNSC_OV_SKCM_BRCA_selectedTumor_filters.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## BRCA noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/filters/filters.Rmd", 
    params = list(
        metadata_tissue = c(BRCA="../../metadata/BRCA/BRCA_technical_metadata.txt"),
        metadata_comp = c(COAD="../../metadata/COAD/COAD_technical_metadata.txt", 
                            GBM="../../metadata/GBM/GBM_technical_metadata.txt", 
                            LUSC="../../metadata/LUSC/LUSC_technical_metadata.txt", 
                            HNSC="../../metadata/HNSC/HNSC_technical_metadata.txt", 
                            OV="../../metadata/OV/OV_technical_metadata.txt",
                            READ="../../metadata/READ/READ_technical_metadata.txt",
                            SKCM="../../metadata/SKCM/SKCM_technical_metadata.txt",
                            BRCA="../../metadata/BRCA/BRCA_technical_metadata.txt"),
        taxa_tissue = c(BRCA="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_comp = c(COAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        GBM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        HNSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        OV="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        READ="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        SKCM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        BRCA="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        feat_tissue = "Project_ID",
        batch_feat = "plate_id",
        output_folder = "../../results/filters/",
        thr_presence = 0.1
    ), 
    output_file = "../../results/filters/BRCA_selectedTumor_vs_COAD_GBM_LUSC_HNSC_OV_READ_SKCM_BRCA_selectedTumor_filters.html"
)

rm(list=ls())
gc(full=TRUE)
