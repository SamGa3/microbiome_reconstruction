
rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## COAD noFFPE CoupledNormalPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

rmarkdown::render("scripts/ml/ml_lasso_classifier.Rmd", 
    params = list(
        metadata = "../../metadata/COAD/COAD_clinical_metadata.txt",
        new_property = "",
        taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedCoupledTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
        match_metadata_to_taxa = "file_id",
        cat_properties = "sample_type",
        values_not_considered = list("unknown"),
        mutate_cat_properties = list(adj_sample_type=list("Primary Tumor"="", "Solid Tissue Normal"="")),
        paired = "patient_id",
        filter = c("sd", 200), 
        total_taxa = "../../data/all_bacteria_species.txt"
    ), 
    output_file = "../../results/ml/COAD_CoupledTumorNormal_ml_lasso.html"
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## COAD noFFPE OnlyPrimary NoRiboZeroG

#---------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/ml/ml_lasso_classifier.Rmd", 
    params = list(
        metadata = "../../metadata/COAD/COAD_clinical_metadata.txt",
        match_metadata_to_taxa = "file_id",
        taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        join = "",
        new_property = "",
        cat_properties = c("MSI_status", "side", "CMS"),
        values_not_considered = list("unknown", "unknown", c("unknown", "NOLBL")),
        mutate_cat_properties = list(adj_MSI=list(high="", low=""), 
                                    adj_side=list(left="", right=""), 
                                    adj_CMS1=list(CMS1="CMS1", other=c("CMS2", "CMS3", "CMS4"))
                                ),
        paired = rep(NULL, 3),
        filter = c("sd", 500),
        total_taxa = "../../data/all_bacteria_species.txt"
    ),
    output_file = "../../results/ml/COAD_selectedTumor_ml_lasso_clinical1.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/ml/ml_lasso_classifier.Rmd", 
    params = list(
        metadata = "../../metadata/COAD/COAD_clinical_metadata.txt",
        match_metadata_to_taxa = "file_id",
        taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        join = "",
        new_property = "",
        cat_properties = rep("CMS", 4),
        values_not_considered = list(c("unknown", "NOLBL"), c("unknown", "NOLBL"), c("unknown", "NOLBL"), c("unknown", "NOLBL")),
        mutate_cat_properties = list(adj_CMS1=list(CMS1="CMS1", other=c("CMS2", "CMS3", "CMS4")), 
                                    adj_CMS2=list(CMS2="CMS2", other=c("CMS1", "CMS3", "CMS4")), 
                                    adj_CMS3=list(CMS3="CMS3", other=c("CMS1", "CMS2", "CMS4")), 
                                    adj_CMS4=list(CMS4="CMS4", other=c("CMS1", "CMS2", "CMS3"))
                                ),
        paired = c(NULL, NULL),
        filter = c("sd", 500),
        total_taxa = "../../data/all_bacteria_species.txt"
    ),
    output_file = "../../results/ml/COAD_selectedTumor_ml_lasso_CMS.html"
)

rm(list=ls())
gc(full=TRUE)
