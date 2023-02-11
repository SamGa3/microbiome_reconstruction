#---------------------------------------------------------------------------------------------------------------------------

## GENE SYMBOL CONVERTER

#---------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/gene_expression/gene_symbol_converter.Rmd", 
  params=list(
    tpm = "../../data/RNAseq/TPM/tcga/TCGA-COAD_tpm.txt",
    column_name = "ensembl_id",
    input_gene = "gene_id",
    output_gene = "gene_name",
    converter_tab = "../../data/RNAseq/TPM/tcga/gene_annotation_v22_gene_length.txt",
    output = c("../../data/RNAseq/TPM/COAD_TPM_tot.txt")
  ), 
  output_file = "../../data/RNAseq/TPM/tcga/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

# names from tiziano, not chosen
rmarkdown::render("scripts/gene_expression/gene_symbol_converter.Rmd", 
  params=list(
    tpm = c("../../data/RNAseq/TPM/COAD_TPM_tot.txt"),
    column_name = "ensembl_id",
    input_gene = "Approved_symbol",
    output_gene = "gene_names",
    converter_tab = "gene_names_from_tiziano_COAD_READ_OnlyTum.txt",
    output = c("../../data/RNAseq/TPM/COAD_TPM_tiziano.txt")
  ), 
  output_file = "../../data/RNAseq/TPM/tcga/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/gene_expression/gene_symbol_converter.Rmd", 
  params=list(
    tpm = "../../data/RNAseq/TPM/tcga/TCGA-GBM_tpm_df.txt",
    column_name = "ensembl_id",
    input_gene = "gene_id",
    output_gene = "gene_name",
    converter_tab = "../../data/RNAseq/TPM/tcga/gene_annotation_v22_gene_length.txt",
    output = c("../../data/RNAseq/TPM/GBM_TPM_tot.txt")
  ), 
  output_file = "../../data/RNAseq/TPM/tcga/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/gene_expression/gene_symbol_converter.Rmd", 
  params=list(
    tpm = "../../data/RNAseq/TPM/tcga/TCGA-LUAD_tpm_df.txt",
    column_name = "ensembl_id",
    input_gene = "gene_id",
    output_gene = "gene_name",
    converter_tab = "../../data/RNAseq/TPM/tcga/gene_annotation_v22_gene_length.txt",
    output = c("../../data/RNAseq/TPM/LUAD_TPM_tot.txt")
  ), 
  output_file = "../../data/RNAseq/TPM/tcga/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/gene_expression/gene_symbol_converter.Rmd", 
  params=list(
    tpm = "../../data/RNAseq/TPM/tcga/TCGA-LUSC_tpm_df.txt",
    column_name = "ensembl_id",
    input_gene = "gene_id",
    output_gene = "gene_name",
    converter_tab = "../../data/RNAseq/TPM/tcga/gene_annotation_v22_gene_length.txt",
    output = c("../../data/RNAseq/TPM/LUSC_TPM_tot.txt")
  ), 
  output_file = "../../data/RNAseq/TPM/tcga/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/gene_expression/gene_symbol_converter.Rmd", 
  params=list(
    tpm = "../../data/RNAseq/TPM/tcga/TCGA-HNSC_tpm_df.txt",
    column_name = "ensembl_id",
    input_gene = "gene_id",
    output_gene = "gene_name",
    converter_tab = "../../data/RNAseq/TPM/tcga/gene_annotation_v22_gene_length.txt",
    output = c("../../data/RNAseq/TPM/HNSC_TPM_tot.txt")
  ), 
  output_file = "../../data/RNAseq/TPM/tcga/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/gene_expression/gene_symbol_converter.Rmd", 
  params=list(
    tpm = "../../data/RNAseq/TPM/tcga/TCGA-OV_tpm_df.txt",
    column_name = "ensembl_id",
    input_gene = "gene_id",
    output_gene = "gene_name",
    converter_tab = "../../data/RNAseq/TPM/tcga/gene_annotation_v22_gene_length.txt",
    output = c("../../data/RNAseq/TPM/OV_TPM_tot.txt")
  ), 
  output_file = "../../data/RNAseq/TPM/tcga/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/gene_expression/gene_symbol_converter.Rmd", 
  params=list(
    tpm = "../../data/RNAseq/TPM/tcga/TCGA-SKCM_tpm_df.txt",
    column_name = "ensembl_id",
    input_gene = "gene_id",
    output_gene = "gene_name",
    converter_tab = "../../data/RNAseq/TPM/tcga/gene_annotation_v22_gene_length.txt",
    output = c("../../data/RNAseq/TPM/SKCM_TPM_tot.txt")
  ), 
  output_file = "../../data/RNAseq/TPM/tcga/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/gene_expression/gene_symbol_converter.Rmd", 
  params=list(
    tpm = "../../data/RNAseq/TPM/tcga/TCGA-READ_tpm_df.txt",
    column_name = "ensembl_id",
    input_gene = "gene_id",
    output_gene = "gene_name",
    converter_tab = "../../data/RNAseq/TPM/tcga/gene_annotation_v22_gene_length.txt",
    output = c("../../data/RNAseq/TPM/READ_TPM_tot.txt")
  ), 
  output_file = "../../data/RNAseq/TPM/tcga/tpm.html"
)

rm(list=ls())
gc(full=TRUE)

rmarkdown::render("scripts/gene_expression/gene_symbol_converter.Rmd", 
  params=list(
    tpm = "../../data/RNAseq/TPM/tcga/TCGA-BRCA_tpm_df.txt",
    column_name = "ensembl_id",
    input_gene = "gene_id",
    output_gene = "gene_name",
    converter_tab = "../../data/RNAseq/TPM/tcga/gene_annotation_v22_gene_length.txt",
    output = c("../../data/RNAseq/TPM/BRCA_TPM_tot.txt")
  ), 
  output_file = "../../data/RNAseq/TPM/tcga/tpm.html"
)
