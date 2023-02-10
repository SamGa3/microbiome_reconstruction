
#---------------------------------------------------------------------------------------------------------------------------

## TCGA

#---------------------------------------------------------------------------------------------------------------------------

## COAD

# noFFPE NormalPrimaryTumor NoRiboZeroG
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name"),
    selection = list(c("Primary Tumor", "Solid Tissue Normal"), "NO", "Illumina TruSeq"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# onlycoupled noFFPE NormalPrimaryTumor NoRiboZeroG
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("coupled", "sample_type", "is_ffpe", "library_name"),
    selection = list(TRUE, c("Primary Tumor", "Solid Tissue Normal"), "NO", "Illumina TruSeq"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedCoupledTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# solid tissue normal
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type"),
    selection = list("Solid Tissue Normal"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# onlycoupled noFFPE OnlyPrimaryTumor NoRiboZeroG
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("coupled", "sample_type", "is_ffpe", "library_name"),
    selection = list(TRUE, "Primary Tumor", "NO", "Illumina TruSeq"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedCoupledTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# noFFPE OnlyPrimaryTumor NoRiboZeroG
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown")),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_norm/COAD/COAD_bacteria_species_merged_unamb_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown")),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_ambig_norm/COAD/COAD_bacteria_species_merged_ambig_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown")),
    output = "../../data/RNAseq/bacteria/raw/merged_ambig_norm/COAD/COAD_selectedTumor_bacteria_species_merged_ambig_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# noFFPE OnlyPrimaryTumor NoRiboZeroG read_length48
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name", "read_length"),
    selection = list("Primary Tumor", "NO", "Illumina TruSeq", "48"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumorReadLength48_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# noFFPE OnlyPrimaryTumor NoRiboZeroG read_length76
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name", "read_length"),
    selection = list("Primary Tumor", "NO", "Illumina TruSeq", "76"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumorReadLength76_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# genus
# noFFPE OnlyPrimaryTumor NoRiboZeroG
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_bacteria_genus_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown")),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_genus_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

## GBM

# OnlyPrimaryTumor NoDupl
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/GBM/GBM_technical_metadata.txt", "../../metadata/GBM/GBM_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "primary_duplicated"),
    selection = list("Primary Tumor", FALSE),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# NormalPrimary NoDupl
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/GBM/GBM_technical_metadata.txt", "../../metadata/GBM/GBM_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "primary_duplicated"),
    selection = list(c("Primary Tumor", "Solid Tissue Normal"), FALSE),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

## LUAD

# noFFPE OnlyPrimaryTumor NoRiboZeroG
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/LUAD/LUAD_technical_metadata.txt", "../../metadata/LUAD/LUAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name"),
    selection = list("Primary Tumor", "NO", "Illumina TruSeq"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# noFFPE NormalPrimary NoRiboZeroG
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/LUAD/LUAD_technical_metadata.txt", "../../metadata/LUAD/LUAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name"),
    selection = list(c("Primary Tumor", "Solid Tissue Normal"), "NO", "Illumina TruSeq"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

## LUSC

# OnlyPrimaryTumor NoDupl
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/LUSC/LUSC_technical_metadata.txt", "../../metadata/LUSC/LUSC_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "primary_duplicated"),
    selection = list("Primary Tumor", FALSE),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# NormalPrimary NoDupl
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/LUSC/LUSC_technical_metadata.txt", "../../metadata/LUSC/LUSC_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "primary_duplicated"),
    selection = list(c("Primary Tumor", "Solid Tissue Normal"), FALSE),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

## HNSC

# noFFPE OnlyPrimaryTumor
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/HNSC/HNSC_technical_metadata.txt", "../../metadata/HNSC/HNSC_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe"),
    selection = list("Primary Tumor", "NO"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# noFFPE NormalPrimary
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/HNSC/HNSC_technical_metadata.txt", "../../metadata/HNSC/HNSC_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe"),
    selection = list(c("Primary Tumor", "Solid Tissue Normal"), "NO"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

## OV

# OnlyPrimaryTumor nomirVana
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/OV/OV_technical_metadata.txt", "../../metadata/OV/OV_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "experimental_protocol_type"),
    selection = list("Primary Tumor", "Allprep RNA Extraction"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

## READ

# OnlyPrimaryTumor
 
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/READ/READ_technical_metadata.txt", "../../metadata/READ/READ_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c(c("sample_type")),
    selection = list("Primary Tumor"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# NormalPrimary
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/READ/READ_technical_metadata.txt", "../../metadata/READ/READ_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c(c("sample_type")),
    selection = list(c("Primary Tumor", "Solid Tissue Normal")),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

## SKCM

# OnlyPrimaryTumor
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/SKCM/SKCM_technical_metadata.txt", "../../metadata/SKCM/SKCM_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type"),
    selection = list("Primary Tumor"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

#  NormalPrimaryTumor
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/SKCM/SKCM_technical_metadata.txt", "../../metadata/SKCM/SKCM_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type"),
    selection = list(c("Primary Tumor", "Solid Tissue Normal")),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

## BRCA

# noFFPE OnlyPrimaryTumor NoRiboZeroG
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/BRCA/BRCA_technical_metadata.txt", "../../metadata/BRCA/BRCA_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name"),
    selection = list("Primary Tumor", "NO", "Illumina TruSeq"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# noFFPE NormalPrimary NoRiboZeroG
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/BRCA/BRCA_technical_metadata.txt", "../../metadata/BRCA/BRCA_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name"),
    selection = list(c("Primary Tumor", "Solid Tissue Normal"), "NO", "Illumina TruSeq"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/BRCA/BRCA_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

#---------------------------------------------------------------------------------------------------------------------------

## IEO

#---------------------------------------------------------------------------------------------------------------------------

# OnlySelected
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params=list(
    metadata = c("../../metadata/IEO/IEO_clinical_metadata.txt", "../../metadata/IEO/IEO_technical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("read_num_selection"),
    selection = list(TRUE),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
  )
)

rm(list=ls())
gc(full=TRUE)

# OnlySelected OnlyTumor
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params = list(
    metadata = c("../../metadata/IEO/IEO_clinical_metadata.txt", "../../metadata/IEO/IEO_technical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("read_num_selection", "sample_type"),
    selection = list(TRUE, "TUM"),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)
