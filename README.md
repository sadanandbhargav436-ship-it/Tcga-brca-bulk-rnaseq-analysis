# Bulk RNA-Seq Analysis of Breast Cancer Using TCGA

This project implements a reproducible analysis pipeline for TCGA breast cancer (TCGA-BRCA) bulk RNA-seq data.

The pipeline includes:

- data acquisition using `TCGAbiolinks`
- preprocessing and normalization of HTSeq raw counts
- differential gene expression analysis comparing primary tumor vs normal breast tissue
- exploratory visualizations and enrichment studies
- survival analysis using TCGA clinical data

## Project structure

- `scripts/install_packages.R` — installs the required R packages
- `scripts/01_download_preprocess.R` — queries and downloads TCGA data, and prepares count and clinical datasets
- `scripts/02_analysis.R` — performs normalization, differential expression, enrichment, and survival analysis
- `data/` — directory for downloaded and preprocessed data files
- `results/` — directory for generated result tables and plots

## Required software

- R (version 4.1 or later recommended)
- An internet connection for downloading TCGA data

## Run the pipeline

Option A (automatic installer — Recommended on Windows):

1. Open PowerShell as Administrator.
2. From the project root run:

```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope Process; .\scripts\install_and_run_R.ps1
```

This script will try to detect R, attempt to download & silently install R (fallbacks to manual install if automatic download fails), run the package installer, then execute the pipeline.

Option B (manual R usage):

From the project root in R, run the pipeline directly:

```r
source("scripts/run_pipeline.R")
```

Or run each script individually:

```r
source("scripts/install_packages.R")
source("scripts/01_download_preprocess.R")
source("scripts/02_analysis.R")
```

## Notes

- The workflow is designed for reproducibility; it creates `data/` and `results/` directories if they do not exist.
- The `gdc_manifest.txt` file is available for reference, but the pipeline downloads TCGA data via `TCGAbiolinks` directly.
- If you already have raw data files, update `scripts/01_download_preprocess.R` to load them from `data/` instead of downloading.

## Output

The pipeline saves the following outputs in `results/`:

- `data/BRCA_SE_raw.rds`
- `data/BRCA_counts_raw.rds`
- `data/BRCA_clinical_raw.rds`
- `data/BRCA_counts.rds`
- `data/BRCA_clinical.rds`
- `BRCA_DEA_results.csv`
- `BRCA_DEA_top20.csv`
- `BRCA_volcano.png`
- `BRCA_heatmap_top50.png`
- `BRCA_survival_summary.txt`
- `BRCA_enrichment_results.csv`
