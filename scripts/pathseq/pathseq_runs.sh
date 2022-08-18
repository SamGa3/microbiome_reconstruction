#!/bin/bash

cd /microbiome_reconstruction

# Run Pathseq on 3 examples (4 cores)
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

gatk PathSeqPipelineSpark \
    --input data/RNAseq/input_bam/example2/ERR2756906.bam \
    --filter-bwa-image data/pathseq_tools/pathseq_tutorial/hg19mini.fasta.img \
    --is-host-aligned true \
    --kmer-file data/pathseq_tools/pathseq_tutorial/hg19mini.hss \
    --microbe-fasta data/pathseq_tools/pathseq_tutorial/e_coli_k12.fasta \
    --microbe-bwa-image data/pathseq_tools/pathseq_tutorial/e_coli_k12.fasta.img \
    --taxonomy-file data/pathseq_tools/pathseq_tutorial/e_coli_k12.db \
    --output data/RNAseq/pathseq_output/example2/bam_out.bam \
    --score-metrics data/RNAseq/pathseq_output/example2/scores.txt \
    --filter-metrics data/RNAseq/pathseq_output/example2/metrics.txt \
    --scores-output data/RNAseq/pathseq_output/example2/score_out.txt \
    --spark-runner LOCAL \
    --spark-master local[4]

gatk PathSeqPipelineSpark \
    --input data/RNAseq/input_bam/example3/ERR2756907.bam \
    --filter-bwa-image data/pathseq_tools/pathseq_tutorial/hg19mini.fasta.img \
    --is-host-aligned true \
    --kmer-file data/pathseq_tools/pathseq_tutorial/hg19mini.hss \
    --microbe-fasta data/pathseq_tools/pathseq_tutorial/e_coli_k12.fasta \
    --microbe-bwa-image data/pathseq_tools/pathseq_tutorial/e_coli_k12.fasta.img \
    --taxonomy-file data/pathseq_tools/pathseq_tutorial/e_coli_k12.db \
    --output data/RNAseq/pathseq_output/example3/bam_out.bam \
    --score-metrics data/RNAseq/pathseq_output/example3/scores.txt \
    --filter-metrics data/RNAseq/pathseq_output/example3/metrics.txt \
    --scores-output data/RNAseq/pathseq_output/example3/score_out.txt \
    --spark-runner LOCAL \
    --spark-master local[4]