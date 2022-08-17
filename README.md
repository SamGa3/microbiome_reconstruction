# Microbiome Reconstruction Workflow

## Reconstruction of the condition- and location-specific colon cancer microbiome from human RNA sequencing data

Association between microbes and cancer has been reported repeatedly; however, it is not clear if molecular tumour properties are connected to specific microbial colonization patterns. This is due mainly to the current technical and analytical strategy limitations to characterise tumour-associated bacteria.
Here, we propose an approach to detect bacterial signals in human RNA sequencing data and associate them with the clinical and molecular properties of the tumours. The method was tested on public datasets from The Cancer Genome Atlas and its accuracy was assessed on a new cohort of colorectal cancer patients.
Our analysis shows that intratumoral microbiome composition is correlated with survival, anatomic location, microsatellite instability, consensus molecular subtype and immune cell infiltration in colon tumours. In particular, we find Faecalibacterium prausnitzii, Coprococcus comes, Bacteroides spp., Fusobacterium spp. and Clostridium spp. to be strongly associated with tumour properties.
In conclusions, we implemented an approach to analyse concurrently clinical and molecular properties of the tumour and the composition of the associated microbiome. Our results may improve patient stratification and pave the path for mechanistic studies on microbiota-tumour crosstalk.

## Requirements

This project was developed on a cluster of 12 computing nodes of 28 cores, 128GB of ram (alignment to microbial genomes, second step) and a local computer of 6 cores, 15.5GB of ram. The alignment to microbial genomes is executed by Pathseq from [GATK](https://gatk.broadinstitute.org/hc/en-us) (version 4.0.10.1). All the scripts are written in R (R version 3.6.1 (2019-07-05) -- "Action of the Toes").

To execute all the steps of this workflow, you can choose between manually downloading all the required packages and tools, pulling our docker container (suggested approach) or creating a conda environment from our yml file (conda doesn't provide the gatk version we used, so the second step of this workflow must be executed by another version of gatk or manually downloading the correct version as described below).

### Clone this repository

The first step is to clone this reporitory:
```bash
git clone https://github.com/SamGa3/microbiome_reconstruction.git
```
and then you must create the folder structure:
```bash
cd /microbiome_reconstruction
chmod u+wrx scripts/make_structure.sh
./scripts/make_structure.sh
```

### Manual installation

To download the gatk Pathseq version (more details in the README.md file inside the downloaded folder):
```bash
cd /microbiome_reconstruction
wget https://github.com/broadinstitute/gatk/releases/download/4.0.10.1/gatk-4.0.10.1.zip
```
The required R packages are listed in the requirements.R file. To download all the required R packages:
```bash
cd microbiome_reconstruction
Rscript requirements.R
```
### Docker

To install Docker on your local computer, follow the instructions described [here](https://docs.docker.com/engine/). Remember that you need root permisions to install and run docker. You can then pull the microbiome_reconstruction container:
```bash
docker pull gaiasamb/microbiome_reconstruction
```
### Conda

To install conda on your local computer, follow the instructions described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). As explained [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), you can create your own conda environment with all the required packages from the microbiome_reconstruction.yml file:
```bash
    # IN PROGRESS
    conda env create -f conda/microbiome_reconstruction.yml
```

## Workflow

Following this tutorial you will test the microbial reads extraction from bam files using a toy example and reproduce the figures of the paper for COAD samples. To obtain the same results for another cancer type, you need to adapt the scripts accordingly to input that cancer type specific data provided or your own data. 

### Script overview

The scripts use Rmarkdown and produce an html file that wraps together the figures. When needed for further analyses, tables are also produced.
The structure of this workflow has a scripts folder that contains all the required scripts, categorised by type of analysis. The results of these scripts are stored in a folder with a similar structure in the results folder. Each analysis requires its own functions, written in the functions.R file in each script subfolder, while a general set of function needed by multiple scripts is listed in general_functions.R script.
All the metadata files and the tables produced by the workflow use a "file_id" column as key to link the information of each sample (for TCGA samples it is the file_id provided by the gdc database). To easily recover the patient's barcode, an extra table linking the file_id to the TCGA barcode has been added.

### 1. Setup

All the scripts below must be run from the microbiome_reconstruction folder and are written with relative paths, so you don't need to modify the paths. Please make sure you are running the codes from microbiome_reconstruction folder.

#### If you have chosen Docker

After pulling the docker container as explained above, you must activate it with the interactive way:
```bash
sudo docker run -it -e DISPLAY -v /your/path/to/microbiome_reconstruction/:/microbiome_reconstruction 7a13d2e31810 /bin/bash
```
Where 7a13d2e31810 is the image id of the repository, to check your own image id you can run:
```bash
sudo docker images
```

#### If you have chosen Conda

After creating your enviroment as explained above, you must activate it:
```bash
source activate microbiome_reconstruction_env
```
Remember that step 2 of tyhe workflow is not supported (the available gatk version in conda is not the one needed here). You can manually download gatk as described above or instal another version of gatk like this:
```bash
conda install -c bioconda gatk
```

### 2. Alignment to microbial genomes - Pathseq

This step was run on the IEO cluster using 16 cores with 50GB of memory, taking advange of Singularity and gatk docker image. Here, as an example, the command to run Pathseq on the file from sample1:
```bash
gatk PathSeqPipelineSpark \
    --input data/RNAseq/input_bam/sample1/input_bam.bam \
    --filter-bwa-image data/pathseq_tools/pathseq_host.fa.img \
    --is-host-aligned true \
    --kmer-file data/pathseq_tools/pathseq_host.bfi \
    --microbe-fasta data/pathseq_tools/pathseq_microbe.fa \
    --microbe-bwa-image data/pathseq_tools/pathseq_microbe.fa.img \
    --taxonomy-file data/pathseq_tools/pathseq_taxonomy.db \
    --output data/RNAseq/pathseq_output/sample1/bam_out.bam \
    --score-metrics data/RNAseq/pathseq_output/sample1/scores.txt \
    --filter-metrics data/RNAseq/pathseq_output/sample1/metrics.txt \
    --scores-output data/RNAseq/pathseq_output/sample1/score_out.txt \
    --spark-runner LOCAL \
    --spark-master local[16]
```
Since this step requires a lot of time, we provide here a toy sample to test Pathseq on your local computer. This example was taken from [gatk Pathseq tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035889911--How-to-Run-the-Pathseq-pipeline). The first step is to download the example references: 
```bash
cd /microbiome_reconstruction/data/pathseq_tools
wget 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/tutorials/datasets/tutorial_10913.tar.gz' -P /microbiome_reconstruction/data/pathseq_tools
tar –xvzf tutorial_10913.tar.gz
```
Three samples were downloaded from [Zmora et al.](https://www.sciencedirect.com/science/article/pii/S0092867418311024?via%3Dihub) from [ENA](https://www.ebi.ac.uk/ena/browser/home) with the ERR2756905-7 accession number and aligned with STAR-2.7.7a following [this](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/) pipeline. We selected the unmapped reads with samtools. The resulting files are located in /microbiome_reconstruction/data/RNAseq/input_bam/ folder. Then you can run the example (few minutes per sample):
```bash
cd /microbiome_reconstruction
gatk PathSeqPipelineSpark \
    --input data/RNAseq/input_bam/example1/ERR2756905.bam \
    --filter-bwa-image data/pathseq_tools/pathseq_tutorial/hg19mini.fasta.img \
    --is-host-aligned true \
    --kmer-file data/pathseq_tools/pathseq_tutorial/hg19mini.hss \
    --microbe-fasta data/pathseq_tools/pathseq_tutorial/e_coli_k12.fasta \
    --microbe-bwa-image data/pathseq_tools/pathseq_tutorial/e_coli_k12.fasta.img \
    --taxonomy-file data/pathseq_tools/pathseq_tutorial/e_coli_k12.db \
    --output data/RNAseq/pathseq_output/example1/bam_out.bam \
    --score-metrics data/RNAseq/pathseq_output/example1/scores.txt \
    --filter-metrics data/RNAseq/pathseq_output/example1/metrics.txt \
    --scores-output data/RNAseq/pathseq_output/example1/score_out.txt \
    --spark-runner LOCAL \
    --spark-master local[4]
```

### 3. Create unambiguous read and score microbiome tables

From this step on, we will run R scripts only, so you need to activate R from the microbiome_reconstruction folder.
```bash
cd /microbiome_reconstruction
../R-3.6.1/bin/R
```

This step of the pipeline creates a table of unambiguous reads and scores of all samples analysed by Pathseq. The results of this step are already available in the final folders (/microbiome_reconstruction/data/RNAseq/bacteria/raw/unamb and /microbiome_reconstruction/data/RNAseq/bacteria/score) for all the analysed TCGA cancer types and IEO cohort samples. Here we apply the pipeline to the toy samples previously analysed to test it on your local computer. 

```R
rmarkdown::render("scripts/microbes_values/microbiome_table_maker.Rmd", 
    params = list(
        sample_sheet = "../../metadata/example/gdc_sample_sheet_fac_simile.tsv",
        input_folder = "../../data/RNAseq/pathseq_output/",
        output_table = "../../data/RNAseq/bacteria/raw/unamb/example_bacteria_species_unamb.txt",
        kingdom_level = "Bacteria",
        taxon_level = "species",
        type_of_value = "unambiguous"
    )
)
rmarkdown::render("scripts/microbes_values/microbiome_table_maker.Rmd", 
    params = list(
        sample_sheet = "../../metadata/example/gdc_sample_sheet_fac_simile.tsv",
        input_folder = "../../data/RNAseq/pathseq_output/",
        output_table = "../../data/RNAseq/bacteria/raw/score/example_bacteria_species_score.txt",
        kingdom_level = "Bacteria",
        taxon_level = "species",
        type_of_value = "score"
    )
)
```

### 4. Create unambiguous score microbiome table

This step of the pipeline creates the table of unambiguous scores. To do this you need both the table of unambiguous reads and score in the data folder. These two tables have been created summarising each sample table output of pathseq (step 3).
As suggested by Pathseq, some taxonomic_id must be merged together. The list of these taxa is in the merged.dmp file, provided with Pathseq.

```R
rmarkdown::render("scripts/microbes_values/microbiome_estimation.Rmd", 
    params = list(
        unamb="../../data/RNAseq/bacteria/raw/unamb/example_bacteria_species_unamb.txt",
        score="../../data/RNAseq/bacteria/raw/score/example_bacteria_species_score.txt",
        unamb_score_norm_tab="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/example_bacteria_species_merged_unamb_score_norm.txt"
    )
)
```
To obtain all the tables needed to reproduce this paper results, you can run this file that wraps up all the runs needed.

```bash
../R-3.6.1/bin/Rscript scripts/microbes_values/microbiome_estimation_commands.R
```

### 5. Selection of samples

This step selects the samples you are interested in given their properties (e.g. only primary tumour, by read length, or to remove duplicates) and creates a new table with the subset of the analysed samples. 
Here we show how to select only primary, not duplicated COAD samples. 

```R
rmarkdown::render("scripts/microbes_values/sample_selection.Rmd", 
  params = list(
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
    match_metadata = "file_id",
    taxa_tab = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_bacteria_species_merged_unamb_score_norm.txt",
    match_taxa = "rownames",
    properties = c("sample_type", "is_ffpe", "library_name"),
    selection = list("Primary Tumor", "NO", c("Illumina TruSeq", "unknown")),
    output = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
  )
)
```
A list of all the subset of samples used in this paper is available in scripts/microbes_values/sample_selection_commands.R and can be run altogether:

```bash
../R-3.6.1/bin/Rscript scripts/microbes_values/sample_selection_commands.R
```

### 6. Technical batch effect

#### Technical batch effect detection

The bacterial signal detected can be affected by technical issues. This step detects the strongest batch effect by selecting the technical property from TCGA metadata that better describes the differences between samples for each cancer type:

```R
rmarkdown::render("scripts/technical_batch_effect/batch_detection.Rmd", 
  params = list(
    tissues = c("COAD", "GBM", "LUAD", "LUSC", "HNSC", "OV", "READ", "SKCM"),
    metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/GBM/GBM_technical_metadata.txt",
                  "../../metadata/LUAD/LUAD_technical_metadata.txt", "../../metadata/LUSC/LUSC_technical_metadata.txt", 
                  "../../metadata/HNSC/HNSC_technical_metadata.txt", "../../metadata/OV/OV_technical_metadata.txt", 
                  "../../metadata/READ/READ_technical_metadata.txt", "../../metadata/SKCM/SKCM_technical_metadata.txt"),
    new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
    taxa = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                  "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
                )
  ), 
  output_file = "../../results/technical_batch_effect/COAD_LUAD_LUSC_HNSC_OV_READ_SKCM_bacteria_species_merged_unamb_score_norm_batch_detection.html"
)
```
A list of all the test done in this paper is available in scripts/technical_batch_effect/batch_effect_detection_commands.R and can be run altogether:

```bash
../R-3.6.1/bin/Rscript scripts/technical_batch_effect/batch_effect_detection_commands.R
```

#### Technical batch effect correction

This step corrects the batch effect detected in the previous step. Here we show how to correct the batch effect due to the plate_id (corrected to sum together plate_ids with too few samples) in COAD:
```R
rmarkdown::render("scripts/technical_batch_effect/batch_correction.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt"),
        join = "columns", 
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        property = "corr_plate_id",
        output = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
    ), 
    output_file = "../../results/technical_batch_effect/COAD_selectedTumor_ComBat_batch_correction_corr_plate_id.html"
)
```
A list of all the correction used in this paper is in scripts/technical_batch_effect/batch_correction_commands.R and can be run altogether:
```bash
../R-3.6.1/bin/Rscript scripts/technical_batch_effect/batch_correction_commands.R
```

#### Technical batch effect comparison

To compare the differences of bacteria estimations before and after the batch correction of plate_id in COAD, LUAD, LUSC, HNSC, OV, READ and SKCM samples:
```R
rmarkdown::render("scripts/technical_batch_effect/batch_correction_comparison.Rmd", 
    params = list(
        tissues = c("COAD", "LUAD", "LUSC", "HNSC", "OV", "READ", "SKCM"),
        metadata = list("../../metadata/COAD/COAD_technical_metadata.txt",
                    "../../metadata/LUAD/LUAD_technical_metadata.txt", 
                    "../../metadata/LUSC/LUSC_technical_metadata.txt", 
                    "../../metadata/HNSC/HNSC_technical_metadata.txt", 
                    "../../metadata/OV/OV_technical_metadata.txt", 
                    "../../metadata/READ/READ_technical_metadata.txt", 
                    "../../metadata/SKCM/SKCM_technical_metadata.txt"),
        taxa_raw = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_corrected = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/LUSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/HNSC_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/OV_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/SKCM_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        batches = rep("plate_id", 7)
    ), 
    output_file = "../../results/technical_batch_effect/COAD_LUAD_LUSC_HNSC_OV_READ_SKCM_ComBat_corr_plate_id_bacteria_species_batch_comparison_plate_id.html"
)
```
Since COAD and READ show a strong batch effect due to the read length used to sequence the samples (48 and 76bp), here we test the differences between samples with different read_length before and after the batch correction of plate_id in COAD and READ samples:
```R
rmarkdown::render("scripts/technical_batch_effect/batch_correction_comparison.Rmd", 
    params = list(
        tissues = c("COAD", "READ"),
        metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/READ/READ_technical_metadata.txt"),
        taxa_raw = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/READ/READ_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_corrected = c("../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                            "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/READ_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        batches = rep("read_length", 2)
    ), 
    output_file = "../../results/technical_batch_effect/COAD_READ_ComBat_corr_plate_id_bacteria_species_batch_comparison_read_length.html"
)
```
A list of all the comparisons done in this paper is in scripts/technical_batch_effect/batch_effect_comparison_commands.R and can be run altogether:
```bash
../R-3.6.1/bin/Rscript scripts/technical_batch_effect/batch_effect_comparison_commands.R
```

### 7. Property association

This represents one of the most important steps of the workflow. It applies PCA to the reconstructed microbiome and associates the PCs to the properties of the tumour. Here the script to find the associations between the COAD reconstructed microbiome and the clinical tumour properties.
```R
rmarkdown::render("scripts/property_association/diversity.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_clinical_metadata.txt"),
        join = "columns",
        taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        cat_properties = c("gender", "bmi", "stage", "CMS", "history_of_other_malignancy", "side", "MSI_status", "CIMP_status", "history_colon_polyps"),
        values_not_considered = list("unknown", "unknown", "unknown", c("unknown", "NOLBL"), c("unknown","inconsistency"), "unknown",
                                       "unknown", "unknown", "unknown"),
        cont_properties = c("percent_normal_cells", "age", "mutation_burden", "stemness", 
                            "aneuploidy_score"),
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        rotations_path="../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_rotations.txt",
        pca_matrix_path = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_pca_tab.txt",
        palette = c("default", "default", "default", "CMS", rep("default", 10))
    ), 
    output_file = "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/COAD_ComBat_plate_id_selectedTumor_property_association.html"
)
```
A list of all the tests used in this paper is in scripts/property_association/porperty_association_commands.R and can be run altogether:
```bash
../R-3.6.1/bin/Rscript scripts/property_association/property_association_commands.R
```

### 8. Survival analysis

This step looks for associations between the reconstructed microbiome and the survival of patients. A list of all the tests used in this paper is in scripts/survival_analysis/survival_analysis_commands.R and can be run altogether:
```bash
../R-3.6.1/bin/Rscript scripts/survival_analysis/survival_analysis_commands.R
```

#### Cox proportional-hazards model

Here we test if any PCs of COAD reconstructed microbiome is associated to the DFS of patients:
```R
rmarkdown::render("scripts/survival_analysis/cox_analysis.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt",
                        "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_pca_tab.txt"),
        join = c("columns"),
        taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        surv_data = "../../metadata/COAD/COAD_cBioPortal_disease_free_survival_firehose.txt",
        survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
        categorical_covariates = "",
        values_not_considered = "",
        timerange_cat = "",
        numeric_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"),
        timerange_cont = list(c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5), c(0, 5))
    ), 
    output_file = "../../results/survival_analysis/COAD_selectedTumor_COX_DFS.html"
)
```

#### Kaplan-Meyer

Since PC4 has been found associated to DFS in Cox proportional-hazard model test and other COAD tumour properties are associated to PC4, we check if PC4 or any og these properties are associated to the DFS of patients:
```R
rmarkdown::render("scripts/survival_analysis/km_analysis.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt", 
                        "../../results/property_association/bacteria_species/ComBat_batch_corrected_plate_id/merged_unamb_score_norm/COAD/tables/COAD_ComBat_corr_plate_id_selectedTumor_pca_tab.txt"),
        join = "columns",
        taxa = "../../data/RNAseq/bacteria/ComBat_plate_id/merged_unamb_score_norm/COAD_ComBat_corr_plate_id_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        surv_data = "../../metadata/COAD/COAD_cBioPortal_disease_free_survival_firehose.txt",
        survival_analysis = c("DFS_YEARS", "patient_status", "DFS"),
        categorical_covariates = c("history_colon_polyps"),
        values_not_considered = c("unknown"),
        timerange_cat = list(c(0, 5)),
        numeric_covariates = c("PC4", "age", "mutation_burden"),
        timerange_cont = list(c(0, 5), c(0, 5), c(0, 5)),
        break_val = 1
    ), 
    output_file = "../../results/survival_analysis/COAD_selectedTumor_KM_DFS_PC4.html"
)
```

### 9. Bacterial species filters

This step selects the species that satisfy the criteria needed to either tissue specificity, not-batch association or most prevalence in tumor COAD reconstructed microbiomes.
```R
rmarkdown::render("scripts/filters/filters.Rmd", 
    params = list(
        metadata_tissue = c(COAD="../../metadata/COAD/COAD_technical_metadata.txt"),
        metadata_comp = c(GBM="../../metadata/GBM/GBM_technical_metadata.txt", 
                            LUAD="../../metadata/LUAD/LUAD_technical_metadata.txt", 
                            LUSC="../../metadata/LUSC/LUSC_technical_metadata.txt", 
                            HNSC="../../metadata/HNSC/HNSC_technical_metadata.txt", 
                            OV="../../metadata/OV/OV_technical_metadata.txt",
                            SKCM="../../metadata/SKCM/SKCM_technical_metadata.txt"),
        taxa_tissue = c(COAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        taxa_comp = c(GBM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/GBM/GBM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUAD="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUAD/LUAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        LUSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/LUSC/LUSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        HNSC="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/HNSC/HNSC_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        OV="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/OV/OV_selectedTumor_bacteria_species_merged_unamb_score_norm.txt", 
                        SKCM="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/SKCM/SKCM_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        feat_tissue = "Project_ID",
        batch_feat = "plate_id",
        output_folder = "../../results/filters/",
        thr_presence = 0.1
    ), 
    output_file = "../../results/filters/COAD_selectedTumor_vs_GBM_LUAD_LUSC_HNSC_OV_SKCM_selectedTumor_filters.html"
)
```

### 10. Identification of species related to tumour properties

This step identifies the species that are related to the tumour properties previously detected as associated to tumour composition of COAD samples.
```R
rmarkdown::render("scripts/identification_related_species/taxa_compositions.Rmd", 
    params = list(
        metadata = c("../../metadata/COAD/COAD_clinical_metadata.txt", "../../metadata/COAD/COAD_technical_metadata.txt", "../../metadata/COAD/COAD_immuneInfiltrationRelative_pbelow05_metadata.txt"),
        join_by = "columns",
        new_property = list(c(old="plate_id", met="corr_plate_id", new_name="corr_plate_id")),
        taxa = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/COAD/COAD_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
        cat_properties = c("CMS", "side", "MSI_status"),
        cat_correction = rep("corr_plate_id", 3),
        cont_properties = c("Mast_cells_activated", "Mast_cells_resting"),
        cont_correction = rep("corr_plate_id", 2),
        values_not_considered = list(c("unknown", "NOLBL"), "unknown", "unknown"),
        taxa_selection = c("../../results/filters/Presence_more0.1samples_COAD.txt",
                            "../../results/filters/HighMeanVsRest_COAD_vs_GBM_LUAD_LUSC_HNSC_OV_SKCM.txt"),
        taxa_selection_approach="intersect",
        palette = c("CMS", "Left_right", "MSI", "jama", "nejm"),
        table_path = "../../results/identification_related_species/tables/COAD_selectedTumor_3filters"
    ), 
    output_file = "../../results/identification_related_species/COAD_selectTumor_3filters_bacteria_species_compositions.html"
)
```

A list of all the analyses is in scripts/identification_related_species/identification_related_species_commands.R and can be run altogether:
```bash
../R-3.6.1/bin/Rscript scripts/identification_related_species/identification_related_species_commands.R
```

### 11. Microbiome classification

This step classifies the sample properties given their reconstructed microbiome:
```R
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
    output_file = "../../results/ml/COAD_CoupledTumorNormal_ml_lasso200.html"
)
```

A list of all the analyses is in scripts/ml/ml_lasso_classifier_commands.R and can be run altogether:
```bash
../R-3.6.1/bin/Rscript scripts/ml/ml_lasso_classifier_commands.R
```

### 12. Microbial pathway analysis

#### General microbial pathway workflow

!NB: this is the only step that must NOT be run in R!

To detect which microbial pathways are present in the analysed samples, we used the tool HUMAnN 3.0, from Huttenhower lab. Before running the tool, we merged together the samples we are interested in (left or right side of the colon or their random subsets). Here we show an example on the already presented toy samples from Zmora et al. 
In the first step we convert the bam files from Pathseq to fastq files with samtools:
```bash
samtools bam2fq data/RNAseq/pathseq_output/example1/bam_out.bam > data/RNAseq/humann_output/example/ERR2756905_out/fastq_out.fastq
samtools bam2fq data/RNAseq/pathseq_output/example2/bam_out.bam > data/RNAseq/humann_output/example/ERR2756906_out/fastq_out.fastq
samtools bam2fq data/RNAseq/pathseq_output/example3/bam_out.bam > data/RNAseq/humann_output/example/ERR2756907_out/fastq_out.fastq
cat data/RNAseq/humann_output/example/ERR2756905_out/fastq_out.fastq data/RNAseq/humann_output/example/ERR2756906_out/fastq_out.fastq data/RNAseq/humann_output/example/ERR2756907_out/fastq_out.fastq > data/RNAseq/humann_output/example/merged_fastq_out.fastq
```
The tool step was run on the IEO cluster using local computer of 8 cores, 25 GB of ram.
Since the reference databases are big (several GBs), here we show an example with the demo databases, as described in HUMAnN tutorial (https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0), (few minutes per sample):
NB: this step requires several minutes to automatically download the required databases.
```bash
humann --input data/RNAseq/humann_output/example/merged_fastq_out.fastq --output data/RNAseq/humann_output/example
```
Given we are using toy samples and demo databases, the output is empty, meaning that no pathways were detected (as expected).
The second step of HUMAnN 3.0 tool involves the split of stratified and unstratified pathways (not stratified by microbes) and the normalisation of the data from RPKs to CPM:
```bash
# left
humann_split_stratified_table --input data/RNAseq/humann_output/left/COAD_selectedTumor_left_pathabundance.tsv --output data/RNAseq/humann_output/left/
humann_renorm_table --input data/RNAseq/humann_output/left/COAD_selectedTumor_left_pathabundance_unstratified.tsv \
                    --output data/RNAseq/humann_output/left/COAD_selectedTumor_left_pathabundance_unstratified_cpm.tsv \
                    --units cpm \
                    --update-snames
# right
humann_split_stratified_table --input data/RNAseq/humann_output/right/COAD_selectedTumor_right_pathabundance.tsv --output data/RNAseq/humann_output/right/
humann_renorm_table --input data/RNAseq/humann_output/right/COAD_selectedTumor_right_pathabundance_unstratified.tsv \
                    --output data/RNAseq/humann_output/right/COAD_selectedTumor_right_pathabundance_unstratified_cpm.tsv \
                    --units cpm \
                    --update-snames
```

#### Bootstrapping

The selection of the subsets of samples to be merged has been done with this script:
```R
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
```
The scripts to obtain both the left and right subsets of samples are in:
```bash
../R-3.6.1/bin/Rscript scripts/pathway_analysis/boothstrapping_commands.R
```
The workflow listed in the general workflow step is applied on each subset of samples: for privacy reasons, the fastq files of the COAD samples are not reported, while the output table of the runs of HUMAnN 3.0 on the bootstrapped sets of samples are in results/pathway_analysis .
The list of the script to obtain unstratified and normalized tables are in:
```bash
./scripts/pathway_analysis/boothstrapped_samples_management_commands.sh
```
To see the distributions of the pathways from the bootstrapped samples, we ran:
```R
rmarkdown::render("scripts/pathway_analysis/pathway_analysis.Rmd", 
    params=list(
        pathways1 = c(left="../../data/RNAseq/humann_output/left/COAD_selectedTumor_left_pathabundance_unstratified_cpm.tsv"),
        pathways2 = c(right="../../data/RNAseq/humann_output/right/COAD_selectedTumor_right_pathabundance_unstratified_cpm.tsv"),
        random_pathways1 = "../../data/RNAseq/humann_output/random_subset_left/",
        base_name1 = "COAD_selectedTumor_left_50random",
        random_pathways2 = "../../data/RNAseq/humann_output/random_subset_right/",
        base_name2 = "COAD_selectedTumor_right_50random",
        output = "../../results/pathway_analysis/COAD_side_table.tsv"
    ), 
    output_file = "../../results/pathway_analysis/COAD_side_bootstrapping.html"
)
```

### 13. Microbiome reconstruction validation

#### IEO cohort comparison: RNA-Seq vs. FISH

Comparison of FISH quantifications of Akkermansia muciniphila and Faecalibacterium prausnitzii with their estimation from human RNA-seq data from this workflow
```R
rmarkdown::render("scripts/comparison/continuous_metadata_analysis.Rmd", 
    params=list(
        metadata1 =  c(FISH="../../metadata/IEO/IEO_FISH_metadata.txt"),
        metadata2 = c(RNAseq="../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"),
        property1 = c("ratio_akk_eub_dapi", "ratio_praus_eub_dapi"),
        property2 = c("239935", "853"),
        taxa_tab = c("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumor_bacteria_species_merged_unamb_score_norm.txt",
                    "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumor_bacteria_species_merged_unamb_score_norm.txt"
                  ),
        output = "../../results/comparison/IEO_correlations_FISH_RNAseq_akk_faec.txt"
    ), 
    output_file = "../../results/comparison/IEO_correlations_FISH_RNAseq_akk_faec.html"
)
```

#### IEO cohort comparison: RNA-Seq vs. 16S

Comparison of 16S quantifications with human RNA-seq estimation from this workflow
```R
rmarkdown::render("scripts/comparison/continuous_metadata_analysis_overview.Rmd", 
    params = list(
        metadata1 =  c("../../metadata/IEO/IEO_16S_metadata.txt"),
        name1 = "16S",
        metadata2 = "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_bacteria_genus_merged_unamb_score_norm.txt",
        name2 = "RNAseq",
        taxa_tab = list("../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt",
                        "../../data/RNAseq/bacteria/raw/merged_unamb_score_norm/IEO/IEO_selectedTumorNormal_bacteria_species_merged_unamb_score_norm.txt"
                    ),
        property1 = c("intersection"),
        property2 = c("intersection"),
        thrs1 = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
        thrs2 = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
        palette="viridis7",
        output = "../../results/comparison/IEO_correlations_16S_RNAseq_"
    ), 
    output_file = "../../results/comparison/IEO_correlations_16S_RNAseq.html"
)
```
