
#---------------------------------------------------------------------------------------------------------------------------

## Create TPM from FPKM sample tables

#---------------------------------------------------------------------------------------------------------------------------

## COAD 

rmarkdown::render("scripts/gene_expression/from_FPKM_to_TPM.Rmd", 
  params=list(
    manifest = "../../data/RNAseq/FPKM/tcga/COAD_fpkm_manifest.tsv",
    dir = "../../data/RNAseq/FPKM/tcga/",
    converter_tab = "../../data/RNAseq/FPKM/tcga/gene_annotation_v22_gene_length.txt",
    output_fpkm = c("../../data/RNAseq/FPKM/"),
    output_tpm = c("../../data/RNAseq/TPM/")
  ), 
  output_file = "../../data/RNAseq/TPM/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

## GBM 

rmarkdown::render("scripts/gene_expression/from_FPKM_to_TPM.Rmd", 
  params=list(
    manifest = "../../data/RNAseq/FPKM/tcga/GBM_fpkm_manifest.tsv",
    dir = "../../data/RNAseq/FPKM/tcga/",
    converter_tab = "../../data/RNAseq/FPKM/tcga/gene_annotation_v22_gene_length.txt",
    output_fpkm = c("../../data/RNAseq/FPKM/"),
    output_tpm = c("../../data/RNAseq/TPM/")
  ), 
  output_file = "../../data/RNAseq/TPM/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

## LUAD 

rmarkdown::render("scripts/gene_expression/from_FPKM_to_TPM.Rmd", 
  params=list(
    manifest = "../../data/RNAseq/FPKM/tcga/LUAD_fpkm_manifest.tsv",
    dir = "../../data/RNAseq/FPKM/tcga/",
    converter_tab = "../../data/RNAseq/FPKM/tcga/gene_annotation_v22_gene_length.txt",
    output_fpkm = c("../../data/RNAseq/FPKM/"),
    output_tpm = c("../../data/RNAseq/TPM/")
  ), 
  output_file = "../../data/RNAseq/TPM/tpm.html"
)

## LUSC 

rmarkdown::render("scripts/gene_expression/from_FPKM_to_TPM.Rmd", 
  params=list(
    manifest = "../../data/RNAseq/FPKM/tcga/LUSC_fpkm_manifest.tsv",
    dir = "../../data/RNAseq/FPKM/tcga/",
    converter_tab = "../../data/RNAseq/FPKM/tcga/gene_annotation_v22_gene_length.txt",
    output_fpkm = c("../../data/RNAseq/FPKM/"),
    output_tpm = c("../../data/RNAseq/TPM/")
  ), 
  output_file = "../../data/RNAseq/TPM/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

## HNSC 

rmarkdown::render("scripts/gene_expression/from_FPKM_to_TPM.Rmd", 
  params=list(
    manifest = "../../data/RNAseq/FPKM/tcga/HNSC_fpkm_manifest.tsv",
    dir = "../../data/RNAseq/FPKM/tcga/",
    converter_tab = "../../data/RNAseq/FPKM/tcga/gene_annotation_v22_gene_length.txt",
    output_fpkm = c("../../data/RNAseq/FPKM/"),
    output_tpm = c("../../data/RNAseq/TPM/")
  ), 
  output_file = "../../data/RNAseq/TPM/tpm.html"
)

## OV 

rmarkdown::render("scripts/gene_expression/from_FPKM_to_TPM.Rmd", 
  params=list(
    manifest = "../../data/RNAseq/FPKM/tcga/OV_fpkm_manifest.tsv",
    dir = "../../data/RNAseq/FPKM/tcga/",
    converter_tab = "../../data/RNAseq/FPKM/tcga/gene_annotation_v22_gene_length.txt",
    output_fpkm = c("../../data/RNAseq/FPKM/"),
    output_tpm = c("../../data/RNAseq/TPM/")
  ), 
  output_file = "../../data/RNAseq/TPM/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

## READ 

rmarkdown::render("scripts/gene_expression/from_FPKM_to_TPM.Rmd", 
  params=list(
    manifest = "../../data/RNAseq/FPKM/tcga/READ_fpkm_manifest.tsv",
    dir = "../../data/RNAseq/FPKM/tcga/",
    converter_tab = "../../data/RNAseq/FPKM/tcga/gene_annotation_v22_gene_length.txt",
    output_fpkm = c("../../data/RNAseq/FPKM/"),
    output_tpm = c("../../data/RNAseq/TPM/")
  ), 
  output_file = "../../data/RNAseq/TPM/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

## SKCM 

rmarkdown::render("scripts/gene_expression/from_FPKM_to_TPM.Rmd", 
  params=list(
    manifest = "../../data/RNAseq/FPKM/tcga/SKCM_fpkm_manifest.tsv",
    dir = "../../data/RNAseq/FPKM/tcga/",
    converter_tab = "../../data/RNAseq/FPKM/tcga/gene_annotation_v22_gene_length.txt",
    output_fpkm = c("../../data/RNAseq/FPKM/"),
    output_tpm = c("../../data/RNAseq/TPM/")
  ), 
  output_file = "../../data/RNAseq/TPM/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

## BRCA 

rmarkdown::render("scripts/gene_expression/from_FPKM_to_TPM.Rmd", 
  params=list(
    manifest = "../../data/RNAseq/FPKM/tcga/BRCA_fpkm_manifest.tsv",
    dir = "../../data/RNAseq/FPKM/tcga/",
    converter_tab = "../../data/RNAseq/FPKM/tcga/gene_annotation_v22_gene_length.txt",
    output_fpkm = c("../../data/RNAseq/FPKM/"),
    output_tpm = c("../../data/RNAseq/TPM/")
  ), 
  output_file = "../../data/RNAseq/TPM/tpm.html"
)
