# Sweet-Corn-Pangenome
Using the Practical Haplotype Graph (PHG) to create a Sweet Corn Pangenome
# 🌽 Sweet Corn PHG Pipeline using PHGv2

This repository provides a full pipeline for building and using a **Practical Haplotype Graph (PHG)** in **sweet corn** using the [PHG version 2](https://github.com/maize-genetics/phg_v2) library developed by Bradbury et al., 2022.

The PHG is a powerful trellis graph representation optimized for crop genomes with high diversity and phased haplotypes. It supports imputation from low-density markers and integrates well with modern breeding tools like **BrAPI**, **rPHG2**, and **rTASSEL**.

---

## 📘 Table of Contents

- [Overview](#overview)
- [Pipeline Steps](#pipeline-steps)
- [Usage](#usage)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Requirements](#requirements)
- [Citation](#citation)
- [Acknowledgments](#acknowledgments)

---

## 📌 Overview

This project builds a **PHG reference for sweet corn**, and uses it to:

- Impute genotypes from short reads via ropeBWT indexing.
- Enable low-cost genomic selection (GS) by imputing from low-density markers.
- Enable GWAS using the imputed genotypes from the PHG path-finding step.

The pipeline follows the official PHGv2 guidelines and includes job scripts for high-performance computing (SLURM).

---

## 🧬 Pipeline Steps

The core PHG workflow includes:

1. **Setup**
   - Initialize conda environment
   - Download and untar the latest PHGv2 release

2. **Reference Construction**
   - `prepare-assemblies` — Clean and format input FASTAs
   - `create-ranges` — Create BED ranges from GFF
   - `align-assemblies` — Align assemblies to the reference
   - `agc-compress` — Compress FASTA sequences
   - `create-ref-vcf` and `create-maf-vcf` — Generate VCF data
   - `load-vcf` — Load VCFs into the PHG TileDB

3. **Imputation**
   - `export-vcf` — Export hVCFs from the DB
   - `rope-bwt-index` — Build ropeBWT index
   - `map-reads` — Map short reads
   - `find-paths` — Impute paths using HMM
   - `hvcf2gvcf` — Convert paths to gVCFs (for downstream use)

---
📂 Input Files
Ensure these are properly formatted and available:

key_file.txt – sample-to-file map

assemblies_list.txt – list of input assemblies

B73v5.fa – sweet corn reference genome (FASTA)

Zm00001eb.1.gff3 – sweet corn genome annotation (GFF)

Raw reads (FASTQ) for imputation

📁 Output Files
Output directories contain:

output/updated_assemblies/ – Renamed FASTAs

output/annotation_B73.bed – BED ranges for PHG

output/B73_alignment/ – MAF alignments

output/vcf_files/ – Reference and MAF VCFs

reads_hvcfs/ – Imputed haplotype VCFs

reads_gvcfs/ – Converted gVCFs

📦 Requirements
PHGv2 ≥ v2.4.x

Conda

Java (≥ 11)

SLURM scheduler

TileDB (via PHG)

📖 Citation
If you use this repository or PHG in your work, please cite:

Bradbury, P. J., et al. (2022). The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation. Bioinformatics. https://doi.org/10.1093/bioinformatics/btac410

🙏 Acknowledgments
This work was developed as part of a PhD project focused on improving genomic selection and GWAS in sweet corn. Special thanks to the PHG development team at Cornell University Maize Genetics Lab. 
