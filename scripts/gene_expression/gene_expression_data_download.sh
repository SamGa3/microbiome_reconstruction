
# COAD
mkdir data/RNAseq/FPKM/tcga/COAD/
../bin/gdc-client download -m data/RNAseq/FPKM/tcga/COAD_fpkm_manifest.tsv -d data/RNAseq/FPKM/tcga/COAD/
gzip -d data/RNAseq/FPKM/tcga/COAD/*/*.gz

# GBM
mkdir data/RNAseq/FPKM/tcga/GBM/
../bin/gdc-client download -m data/RNAseq/FPKM/tcga/GBM_fpkm_manifest.tsv -d data/RNAseq/FPKM/tcga/GBM/
gzip -d data/RNAseq/FPKM/tcga/GBM/*/*.gz

# LUAD
mkdir data/RNAseq/FPKM/tcga/LUAD/
../bin/gdc-client download -m data/RNAseq/FPKM/tcga/LUAD_fpkm_manifest.tsv -d data/RNAseq/FPKM/tcga/LUAD/
gzip -d data/RNAseq/FPKM/tcga/LUAD/*/*.gz

# LUSC
mkdir data/RNAseq/FPKM/tcga/LUSC/
../bin/gdc-client download -m data/RNAseq/FPKM/tcga/LUSC_fpkm_manifest.tsv -d data/RNAseq/FPKM/tcga/LUSC/
gzip -d data/RNAseq/FPKM/tcga/LUSC/*/*.gz

# HNSC
mkdir data/RNAseq/FPKM/tcga/HNSC/
../bin/gdc-client download -m data/RNAseq/FPKM/tcga/HNSC_fpkm_manifest.tsv -d data/RNAseq/FPKM/tcga/HNSC/
gzip -d data/RNAseq/FPKM/tcga/HNSC/*/*.gz

# OV
mkdir data/RNAseq/FPKM/tcga/OV/
../bin/gdc-client download -m data/RNAseq/FPKM/tcga/OV_fpkm_manifest.tsv -d data/RNAseq/FPKM/tcga/OV/
gzip -d data/RNAseq/FPKM/tcga/OV/*/*.gz

# READ
mkdir data/RNAseq/FPKM/tcga/READ/
../bin/gdc-client download -m data/RNAseq/FPKM/tcga/READ_fpkm_manifest.tsv -d data/RNAseq/FPKM/tcga/READ/
gzip -d data/RNAseq/FPKM/tcga/READ/*/*.gz

# SKCM
mkdir data/RNAseq/FPKM/tcga/SKCM/
../bin/gdc-client download -m data/RNAseq/FPKM/tcga/SKCM_fpkm_manifest.tsv -d data/RNAseq/FPKM/tcga/SKCM/
gzip -d data/RNAseq/FPKM/tcga/SKCM/*/*.gz

# BRCA
mkdir data/RNAseq/FPKM/tcga/BRCA/
../bin/gdc-client download -m data/RNAseq/FPKM/tcga/BRCA_fpkm_manifest.tsv -d data/RNAseq/FPKM/tcga/BRCA/
gzip -d data/RNAseq/FPKM/tcga/BRCA/*/*.gz
