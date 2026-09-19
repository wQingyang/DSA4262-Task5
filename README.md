# DSA4262 Task 5 – Long-read RNA-seq Workflow

This repository contains the Nextflow workflow developed for Task 5 of the DSA4262 genomics assignment.

The workflow processes Oxford Nanopore long-read RNA-seq data from the Singapore Nanopore Expression Project (SG-NEx).

## Workflow overview

The workflow performs the following steps:

1. Read alignment using Minimap2
2. SAM-to-BAM conversion, sorting and indexing using Samtools
3. BAM-level quality control
4. Transcript discovery and quantification using Bambu

The workflow was adapted from the long-read RNA-seq Nextflow workflow introduced in DSA4262 Genomics Workshop 3 and extended for the requirements of Task 5.

The main extensions include:

- processing all four samples used in Task 3;
- automatically using protocol-specific Minimap2 parameters for direct RNA and cDNA reads;
- BAM indexing;
- alignment-based quality control;
- support for Bambu analysis with and without genome annotations;
- use of Nextflow caching and `-resume`.

## Input files

The workflow requires:

- Oxford Nanopore FASTQ files;
- GRCh38 reference genome FASTA;
- Ensembl GRCh38 release 91 GTF annotation.

The sequencing and reference data are not included in this repository because of their large file sizes.

## Minimap2 alignment

Different Minimap2 parameters are applied according to the sequencing protocol.

Direct RNA:

```bash
minimap2 -ax splice -uf -k14 reference.fa reads.fastq.gz
```

cDNA:

```bash
minimap2 -ax splice reference.fa reads.fastq.gz
```

## Quality control

QC is performed after alignment using `samtools flagstat`, `samtools stats`, and additional alignment summary statistics.

The workflow records:

- number of primary reads;
- number of mapped primary reads;
- primary-read mapping percentage;
- average read length.

For this assignment, the user-defined QC criteria were:

```text
primary reads >= 1000
mapping rate >= 70%
```

A sample is labelled `PASS` only if both criteria are satisfied.

## Bambu scenarios

Two Bambu analyses are implemented.

### Scenario 1: annotation-guided

Bambu is supplied with the Ensembl GRCh38 release 91 annotations.

### Scenario 2: without genome annotations

Bambu is run without the annotation object. Nextflow is executed using the `-resume` flag so that unchanged alignment, BAM processing and QC steps can be retrieved from the Nextflow cache.

## Output

The workflow generates:

- sorted and indexed BAM files;
- QC summary files;
- Bambu transcript annotations in GTF format;
- transcript read-count tables;
- gene read-count tables;
- serialized Bambu output;
- Nextflow execution reports.

## Task 5 execution

Commands used to execute both workflow scenarios and inspect QC/runtime results are provided in:

```text
run_commands.sh
```

The main workflow is:

```text
task5_longread_rnaseq.nf
```

## Software

The workflow uses:

- Nextflow
- Minimap2
- Samtools
- R
- Bambu

## Notes

Large sequencing data, BAM files, reference genome files, Nextflow work directories and workflow output files are intentionally excluded from this repository.
