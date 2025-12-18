# Cocaine Multiome Analysis Pipeline (Drosophila)

## About the Project

This repository contains shell, R, and Python scripts, along with **Cell Ranger ARC** workflows, used to process and analyze **single-nucleus multiome (RNA + ATAC)** data generated from a six–genetic-background experiment in *Drosophila melanogaster* brain.

The experiment examined **cocaine vs. sucrose exposure**, across **male and female flies**, with the goal of identifying:

- Molecular and signaling pathways responsive to cocaine exposure
- Gene regulatory networks affected by treatment
- Genetic-background– and sex-specific responses

---

## Built With

<p align="left">
  <a href="https://www.10xgenomics.com/support/software/cell-ranger-arc" target="_blank">
    <img src="https://img.shields.io/badge/Cell%20Ranger%20ARC-10xGenomics-blue" alt="Cell Ranger ARC" />
  </a>
  <a href="https://www.r-project.org/" target="_blank">
    <img src="https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white" alt="R" />
  </a>
  <a href="https://www.python.org/" target="_blank">
    <img src="https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white" alt="Python" />
  </a>
  <a href="https://www.gnu.org/software/bash/" target="_blank">
    <img src="https://img.shields.io/badge/Bash-4EAA25?logo=gnu-bash&logoColor=white" alt="Bash" />
  </a>
  <a href="https://snakemake.readthedocs.io/en/stable/" target="_blank">
    <img src="https://img.shields.io/badge/Snakemake-6.4.1-orange" alt="Snakemake" />
  </a>
</p>

---

## Cell Ranger ARC Pipeline

### FASTQ Generation (`mkfastq`)

- **Tool:** `cellranger-arc mkfastq`
- **Purpose:** Convert BCL sequencing data into FASTQ files
- **Location:** `cellranger/mkfastq/`
- **Contents:**
  - `mkfastq_submitter.sh`
  - Sample sheet `.csv` files

These scripts submit `mkfastq` jobs to the cluster to generate FASTQs for downstream processing.

---

### Multiome Quantification (`count`)

- **Tool:** `cellranger-arc count`
- **Workflow:** Snakemake-based
- **Snakemake version:** `6.4.1`

Each batch directory contains:
- `Snakefile`
- YAML configuration file
- SLURM submitter script

Additional requirements:
- A cluster-specific SLURM configuration YAML file
- A `slurm/` directory for logs and job metadata

The individual count summary `.csv` files were concatenated to generate a **project-level summary table**.

---

### Sample Aggregation (`aggr`)

- **Tool:** `cellranger-arc aggr`
- **Purpose:** Aggregate count matrices from all samples
- **Samples:** 48 primary samples (+1 extra replicate)
- **Location:** `cellranger/aggr/`

Contains:
- `aggr_submitter.sh`
- Aggregation `.csv` file

---

## Gene Expression Analysis (GEX)

The `unimodal/` directory reflects analyses where **RNA and ATAC modalities are treated separately**.

### Pseudobulk Differential Expression

**Location:** `unimodal/gex/`

#### Core Script
- `pseudobulk_anova_seurat_avglognorm.R`

This script performs:

- Pseudobulking by cluster
- Log2 normalization
- Quality control filtering
- Metadata factor construction for multifactorial models
- Omnibus ANOVA:
  - Full model: **Line × Treatment × Sex (L×T×S)**
  - Reduced model: **Line × Treatment (L×T)** (sex-averaged)
- Multiple-testing correction of term-level p-values
- Filtering for treatment-associated FDR-significant genes
- Export of differential expression results
- Generation of **LSMeans** per gene, sample, and cluster

**Input:** Cluster number (passed as a command-line argument)

---

#### Supporting Scripts

- `cluster_to_excel.py`
  - Aggregates LSMeans `.txt` files
  - Exports formatted Excel workbooks

- `summarize_table.R`
  - Aggregates and numerically summarizes differential expression results
  - Operates on `_seurat_avglogCPM_filter.txt` files

- `dot_matrix.R`, `dot_matrix_log10.R`
  - Aggregate FDR values across clusters for visualization

- `GEX_only_analysis_Seuratv5.R`
  - Exploratory, unimodal GEX-only analysis
  - Used for QC and cell-type inspection
  - **Not used for formal inference**

---

#### Job Submission

- `sbatch_pseudobulk.sh`
  - SLURM submitter for `pseudobulk_anova_seurat_avglognorm.R`
  - Accepts cluster number as input

This script was executed via a bash loop to process **44 clusters** (48 total, 4 removed for low quality).

---

## ATAC-seq Analysis

The `atac/` directory contains **annotation** and **formal analysis** workflows.

### Peak Annotation

**Location:** `atac/annotation/`

- Peaks were annotated using **GRanges `subsetByOverlaps()`**
- Chosen over TSS-only annotation due to its greater generality

---

### Differential Open Region (DOR) Analysis

**Location:** `atac/analysis/`

#### Core Script
- `pseudobulk_anova_seurat_avglognorm.R`

This script performs:

- Pseudobulking of ATAC counts
- Log2 normalization
- QC filtering
- Multifactorial omnibus ANOVA with **TSS enrichment as a covariate**
- Multiple-testing correction
- Filtering for treatment-associated FDR-significant regions
- Export of DOR results
- LSMeans generation per region and sample
  - Full model: **L×T×S**
  - Reduced model: **L×T**

**Input:** Cluster number

---

#### Supporting Scripts

- `summarize_table.R`
  - Aggregates DOR results from `_seurat_avglogCPM_filter.txt`

- `cluster_to_excel.py`
  - Aggregates LSMeans results into Excel files

---

#### Job Submission Scripts

- `sbatch_pseudobulk.sh`
- `cs_filter_dor.sh`
- `lsmeans.sh`

Each script submits a cluster-specific job to SLURM and routes the cluster number to the corresponding R script.

---

### Peak Comparison and Filtering

**Location:** `atac/analysis/peak_comparison_filtering/`

Purpose:
- Restrict DORs to **reliably detected, cluster-specific peaks**
- Prevent inflation of DOR counts from marginal or spurious signal

Scripts:

- `comp.R`
  - Filters DOR peaks
  - Regenerates per-cluster DOR result files

- `comp_lsmeans_all_peaks.R`
  - Generates LSMeans for **all cluster-specific peaks**
  - Used for data submission

---

## Contact

**Vijay Shankar**  
🔗 https://scienceweb.clemson.edu/ihg/dr-vijay-shankar-2/  
📧 vshanka@clemson.edu

**Project Repository:**  
https://github.com/cfg824157947/Cocaine/

---

## Acknowledgments

- [Clemson Institute for Human Genetics (IHG)](https://scienceweb.clemson.edu/ihg/)
- [COBRE in Human Genetics – P20 GM139769](https://reporter.nih.gov/project-details/11015918)
- [Secretariat HPC](https://secretariat.readthedocs.io/en/latest/)

