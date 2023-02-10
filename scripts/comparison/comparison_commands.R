#---------------------------------------------------------------------------------------------------------------------------

## FISH comaprison to RNA-seq values

#---------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
gc(full=TRUE)

# OnlyTumor, final print pictures
rmarkdown::render("../../microbes/continuous_metadata_analysis.Rmd", params=list(
  metadata1 =  c(FISH="../../metadata/IEO/metadata_FISH.txt"),
  metadata2 = c(RNAseq="../../bacteria/RNAseq/raw/merged_unamb_score_norm/IEO/IEOlist1_noNETOnlySelOnlyTum_bacteria_species_merged_unamb_score_norm.txt"),
  taxa_tab = list(c("../../bacteria/RNAseq/raw/merged_unamb_score_norm/IEO/IEOlist1_noNETOnlySelOnlyTum_bacteria_species_merged_unamb_score_norm.txt"),
                  c("../../bacteria/RNAseq/raw/merged_unamb_score_norm/IEO/IEOlist1_noNETOnlySelOnlyTum_bacteria_species_merged_unamb_score_norm.txt")
                ),
  feature1 = c("mean_ratio_akk_eub_dapi", "mean_ratio_praus_eub_dapi"),
  feature2 = c("239935", "853"),
  group_plots = list(list(metadata1=c(akk="mean_ratio_akk_eub_dapi", praus="mean_ratio_praus_eub_dapi"),
                          metadata2=c(akk="239935", praus="853"),
                          new_labs=c("FISH", "RNAseq", "akk_praus_eub_dapi")
                      )
                  ),
  picture_format = "pdf",
  picture_feat = c("general", "akk_praus_eub_dapi"),
  picture_path = "../../results/comparison/FISH/images/IEO_noNETOnlySelOnlyTumorlist1_correlations_FISH_RNAseq"
), output_file = "../../results/comparison/FISH/IEO_noNETOnlySelOnlyTumorlist1_correlations_FISH_RNAseq_akk_praus.html")

rm(list=ls())
gc(full=TRUE)
