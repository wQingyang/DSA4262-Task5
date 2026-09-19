#!/usr/bin/env bash

# ============================================================
# DSA4262 Task 5
# Long-read RNA-seq workflow execution commands
# ============================================================


# Activate the environment containing Bambu, minimap2,
# samtools, R and the required packages.
conda activate bambu_env


# Move to the Task 3 working directory containing the
# FASTQ and reference files.
cd ~/task3


# ============================================================
# Scenario 1
# Bambu with genome annotations
# ============================================================

nextflow run nextflow/task5_longread_rnaseq.nf \
    --reads "$PWD/fastq/*.fastq.gz" \
    --refFa "$PWD/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa" \
    --refGtf "$PWD/reference/Homo_sapiens.GRCh38.91.gtf" \
    --outdir "$PWD/results" \
    --scenario annotated \
    --with_annotations true \
    -with-report "$PWD/reports/report_annotated.html"


# ============================================================
# Scenario 2
# Bambu without genome annotations
#
# -resume allows unchanged upstream processes to be
# retrieved from the Nextflow cache.
# ============================================================

nextflow run nextflow/task5_longread_rnaseq.nf \
    --reads "$PWD/fastq/*.fastq.gz" \
    --refFa "$PWD/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa" \
    --refGtf "$PWD/reference/Homo_sapiens.GRCh38.91.gtf" \
    --outdir "$PWD/results" \
    --scenario no_annotation \
    --with_annotations false \
    -resume \
    -with-report "$PWD/reports/report_no_annotation.html"


# ============================================================
# Inspect QC summaries
# ============================================================

cat results/qc/*_qc.tsv


# ============================================================
# Inspect Nextflow execution history / runtime
# ============================================================

nextflow log
