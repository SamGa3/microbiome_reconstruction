
#---------------------------------------------------------------------------------------------------------------------------

## select samples

#---------------------------------------------------------------------------------------------------------------------------

## COAD noFFPE OnlyPrimaryTumor NoRiboZeroG

rm(list = ls())
gc(full=TRUE)

# side left
rmarkdown::render("scripts/pathway_analysis/sample_bootstrapping.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_RNAseq_barcodes.txt", "../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"), 
    match_metadata = "file_id",
    properties = c("sample_type", "is_ffpe", "library_name", "side"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown"), "left"),
    column_selected = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    subset_size = 50,
    random_tries = 30,
    output = "../../results/pathway_analysis/bootstrapped_samples/COAD_selectedTumor_left"
  )
)

rm(list = ls())
gc(full=TRUE)

# side right
rmarkdown::render("scripts/pathway_analysis/sample_bootstrapping.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_RNAseq_barcodes.txt", "../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"), 
    match_metadata = "file_id",
    properties = c("sample_type", "is_ffpe", "library_name", "side"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown"), "right"),
    column_selected = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    subset_size = 50,
    random_tries = 30,
    output = "../../results/pathway_analysis/bootstrapped_samples/COAD_selectedTumor_right"
  )
)

rm(list = ls())
gc(full=TRUE)
