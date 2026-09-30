# Drosophila Cocaine Multiome Analysis

Research workflow for analyzing single-cell RNA and ATAC data in a cocaine-related study. The code preprocesses 10x multiome data, recalls cell-type-specific ATAC peaks, constructs LINGER gene-regulatory networks, analyzes differential transcription-factor activity (DTFA), maps GWAS variants to regulatory networks, and creates downstream figures and enrichment results.

## Workflow overview

```text
10x RNA + ATAC matrix
        |
        +-- QC, normalization, clustering, and RNA/ATAC integration
        |
        +-- Cell-type fragment splitting -> MACS2 peak calling -> peak counts
        |
        +-- LINGER GRN -> TF-RE, RE-TG, and TF-TG networks
        |
        +-- DTFA, SNP-RE-TF-TG networks, motif logos, and enrichment analysis
```

Here, **TF** means transcription factor, **RE** means regulatory element, and **TG** means target gene. The analyses use the `dm6` fly genome assembly.

## Repository contents

| Path | Purpose |
| --- | --- |
| `Shscript/Split.sh` | Splits an ATAC fragment table into one BED file per cell type using barcode-to-label mappings. |
| `Shscript/umacs2.sh` | Calls cell-type-specific ATAC peaks with MACS2. |
| `Shscript/calculate_peak_counts.sh` | Merges called peaks and calculates fragment counts over the union peak set with BEDTools. |
| `Shscript/RECALL_PEAK.sh` | SLURM wrapper for the peak-recall workflow. |
| `org_files/preprose.org` | Performs RNA/ATAC quality control, normalization, Harmony integration, clustering, and joint embedding. |
| `org_files/Generate_LINGER_network.org` | Generate LINGER regulatory networks. |
| `org_files/DTFA_and_network_analysis.org` | Calculates differential TF activity and produces network statistics, plots, and top-edge tables. |
| `org_files/SNP_MAP.org` | Converts GWAS variants to BED format and prepares TF-RE-TG networks and regulatory-element BED files. |
| `org_files/RE_SNP_net.org` | Combines SNP-RE overlaps with TF-RE-TG networks and exports Cytoscape node/edge tables. |
| `org_files/Motif_figure.org` | Reads HOMER motif matrices and generates individual or batched sequence-logo figures. |
| `org_files/DIOPT.org` | Filters DIOPT fly-to-human ortholog mappings and performs pathway/enrichment analysis. |

## Requirements

The workflow was designed for a Linux HPC environment. Its main dependencies are:

- Bash, AWK, and standard Unix command-line tools
- SLURM (only for the supplied batch-job directives)
- [MACS2](https://github.com/macs3-project/MACS)
- [BEDTools](https://bedtools.readthedocs.io/)
- Python 3 and Jupyter
- Python packages: `anndata`, `gseapy`, `harmonypy`, `igraph`, `logomaker`, `matplotlib`, `muon`, `numpy`, `pandas`, `psutil`, `scanpy`, `scikit-learn`, `scipy`, `seaborn`, `statsmodels`, `tabulate`, and `torch`
- `LingerGRN` and its required reference data
- Emacs Org mode with Jupyter/Babel support if the `.org` files are run directly

Install `LingerGRN` separately according to the version used for your analysis(recommendation).

## Input data

The repository does not include the study data or external reference files. The workflow expects:

- A 10x HDF5 multiome matrix containing `Gene Expression` and `Peaks` features
- A tab-separated cell-label file with barcode in column 1 and cell type/cluster in column 2
- A tab-separated ATAC fragment file with cell barcode in column 4
- Sample and cell-type metadata used by the DTFA notebooks
- LINGER reference data for `dm6`
- DIOPT ortholog tables and HOMER-format motif matrices for downstream analyses

## Running the workflow

### 1. Configure paths

This is research code and currently contains absolute HPC paths. Before execution, replace the input, output, environment, font, and reference-data paths in both `Shscript/*.sh` and `org_files/*.org` with paths available on your system.

The scripts should be reviewed as templates rather than treated as a portable, end-to-end command. In particular, verify the SLURM partition, memory limits, email address, conda activation path, genome setting, and output directories.

### 2. Recall ATAC peaks

The individual shell stages can be run with:

```bash
bash Shscript/Split.sh condition_label.tsv atac_fragments.tsv cluster_fragments/
bash Shscript/umacs2.sh cluster_fragments/ cluster_peaks/
bash Shscript/calculate_peak_counts.sh
```

Create the output directories first and include a trailing `/` in directory arguments because the current scripts concatenate directory and file names directly. `calculate_peak_counts.sh` uses paths defined inside the script, so edit those variables before running it.

On a configured SLURM cluster, the wrapper can instead be submitted after its paths have been updated:

```bash
sbatch Shscript/RECALL_PEAK.sh
```

### 3. Run the analysis notebooks

Open the Org files in Emacs with a working Jupyter Python kernel (python code blocks)
Several code blocks depend on variables created by earlier blocks in the same Jupyter session.

1. `org_files/preprose.org`
2. `org_files/Generate_LINGER_network.org`
3. `org_files/DTFA_and_network_analysis.org`
4. `org_files/SNP_MAP.org`
5. `org_files/RE_SNP_net.org`
6. `org_files/Motif_figure.org`
7. `org_files/DIOPT.org`


## Main outputs

Depending on the selected notebook blocks, the workflow produces:

- Filtered RNA and ATAC `AnnData` objects
- Target-gene and regulatory-element pseudobulk matrices
- Cell-population TF-RE, RE-TG, and TF-TG network tables
- DTFA statistics, heatmaps, scatterplots, and ranked network edges
- Cell-type-specific Cytoscape node and edge tables
- SNP-to-regulatory-element BED intersections
- Motif-logo figures
- Fly-to-human ortholog enrichment results



## Citation

If you use this workflow, cite the underlying tools and data sources used in your analysis, including Scanpy, Harmony, MACS2, BEDTools, LINGER, DIOPT, HOMER.
