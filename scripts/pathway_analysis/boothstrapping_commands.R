
#---------------------------------------------------------------------------------------------------------------------------

## select samples

#---------------------------------------------------------------------------------------------------------------------------

## COAD noFFPE OnlyPrimaryTumor NoRiboZeroG

rm(list = ls())
gc(full=TRUE)

# side 

# left
rmarkdown::render("scripts/pathway_analysis/sample_bootstrapping.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_RNAseq_barcodes.txt", "../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"), 
    match_metadata = "file_id",
    properties = c("sample_type", "is_ffpe", "library_name", "side"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown"), "left"),
    column_selected = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    subset_size = 55,
    random_tries = 50,
    output = "../../results/pathway_analysis/bootstrapped_samples/COAD_selectedTumor_left"
  )
)

rm(list = ls())
gc(full=TRUE)

# right
rmarkdown::render("scripts/pathway_analysis/sample_bootstrapping.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_RNAseq_barcodes.txt", "../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"), 
    match_metadata = "file_id",
    properties = c("sample_type", "is_ffpe", "library_name", "side"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown"), "right"),
    column_selected = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    subset_size = 72,
    random_tries = 50,
    output = "../../results/pathway_analysis/bootstrapped_samples/COAD_selectedTumor_right"
  )
)

rm(list = ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

# CMS 

# CMS1
rmarkdown::render("scripts/pathway_analysis/sample_bootstrapping.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_RNAseq_barcodes.txt", "../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"), 
    match_metadata = "file_id",
    properties = c("sample_type", "is_ffpe", "library_name", "CMS"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown"), "CMS1"),
    column_selected = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    subset_size = 17,
    random_tries = 50,
    output = "../../results/pathway_analysis/bootstrapped_samples/COAD_selectedTumor_CMS1"
  )
)

rm(list = ls())
gc(full=TRUE)

# CMS234
rmarkdown::render("scripts/pathway_analysis/sample_bootstrapping.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_RNAseq_barcodes.txt", "../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"), 
    match_metadata = "file_id",
    properties = c("sample_type", "is_ffpe", "library_name", "CMS"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown"), c("CMS2", "CMS3", "CMS4")),
    column_selected = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    subset_size = 86,
    random_tries = 50,
    output = "../../results/pathway_analysis/bootstrapped_samples/COAD_selectedTumor_CMS234"
  )
)

rm(list = ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## mutation_burden

# high
rmarkdown::render("scripts/pathway_analysis/sample_bootstrapping.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_RNAseq_barcodes.txt", "../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"), 
    match_metadata = "file_id",
    properties = c("sample_type", "is_ffpe", "library_name", "mutation_burden"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown"), "high"),
    column_selected = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    subset_size = 17,
    random_tries = 50,
    output = "../../results/pathway_analysis/bootstrapped_samples/COAD_selectedTumor_mutation_burden_high"
  )
)

rm(list = ls())
gc(full=TRUE)

# low
rmarkdown::render("scripts/pathway_analysis/sample_bootstrapping.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_RNAseq_barcodes.txt", "../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"), 
    match_metadata = "file_id",
    properties = c("sample_type", "is_ffpe", "library_name", "mutation_burden"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown"), "low"),
    column_selected = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    subset_size = 95,
    random_tries = 50,
    output = "../../results/pathway_analysis/bootstrapped_samples/COAD_selectedTumor_mutation_burden_low"
  )
)

rm(list = ls())
gc(full=TRUE)
