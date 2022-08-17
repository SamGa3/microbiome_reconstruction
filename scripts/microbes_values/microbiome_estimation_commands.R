
#---------------------------------------------------------------------------------------------------------------------------

## TCGA

#---------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
gc(full=TRUE)

# COAD

rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
  params = list(
    ambig = "../../data/RNAseq/bacteria/raw/ambig/COAD_bacteria_species_ambig.txt",
    unamb = "../../data/RNAseq/bacteria/raw/unamb/COAD_bacteria_species_unamb.txt",
    score = "../../data/RNAseq/bacteria/raw/score/COAD_bacteria_species_score.txt",
    ambig_norm_tab = "../../data/RNAseq/bacteria/raw/merged_ambig_norm/COAD/COAD_bacteria_species_merged_ambig_norm.txt",  
    unamb_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_norm/COAD/COAD_bacteria_species_merged_unamb_norm.txt",  
    score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_score_norm/COAD/COAD_bacteria_species_merged_score_norm.txt",  
    unamb_score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# GBM

rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
  params = list(
    unamb = "../../data/RNAseq/bacteria/raw/unamb/GBM_bacteria_species_unamb.txt",
    score = "../../data/RNAseq/bacteria/raw/score/GBM_bacteria_species_score.txt",
    unamb_score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# LUAD

rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
  params = list(
    unamb = "../../data/RNAseq/bacteria/raw/unamb/LUAD_bacteria_species_unamb.txt",
    score = "../../data/RNAseq/bacteria/raw/score/LUAD_bacteria_species_score.txt",
    unamb_score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# LUSC

rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
  params = list(
    unamb = "../../data/RNAseq/bacteria/raw/unamb/LUSC_bacteria_species_unamb.txt",
    score = "../../data/RNAseq/bacteria/raw/score/LUSC_bacteria_species_score.txt",
    unamb_score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# HNSC

rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
  params = list(
    unamb = "../../data/RNAseq/bacteria/raw/unamb/HNSC_bacteria_species_unamb.txt",
    score = "../../data/RNAseq/bacteria/raw/score/HNSC_bacteria_species_score.txt",
    unamb_score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# OV

rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
  params = list(
    unamb = "../../data/RNAseq/bacteria/raw/unamb/OV_bacteria_species_unamb.txt",
    score = "../../data/RNAseq/bacteria/raw/score/OV_bacteria_species_score.txt",
    unamb_score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# READ

rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
  params = list(
    unamb = "../../data/RNAseq/bacteria/raw/unamb/READ_bacteria_species_unamb.txt",
    score = "../../data/RNAseq/bacteria/raw/score/READ_bacteria_species_score.txt",
    unamb_score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# SKCM

rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
  params = list(
    unamb = "../../data/RNAseq/bacteria/raw/unamb/SKCM_bacteria_species_unamb.txt",
    score = "../../data/RNAseq/bacteria/raw/score/SKCM_bacteria_species_score.txt",
    unamb_score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## IEO

#---------------------------------------------------------------------------------------------------------------------------

# IEO

rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
  params = list(
    unamb = "../../data/RNAseq/bacteria/raw/unamb/IEO_bacteria_species_unamb.txt",
    score = "../../data/RNAseq/bacteria/raw/score/IEO_bacteria_species_score.txt",
    unamb_score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
  params = list(
    unamb = "../../data/RNAseq/bacteria/raw/unamb/IEO_bacteria_genus_unamb.txt",
    score = "../../data/RNAseq/bacteria/raw/score/IEO_bacteria_genus_score.txt",
    unamb_score_norm_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_bacteria_genus_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)
