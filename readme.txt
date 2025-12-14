==============================================================================
PD Microglia snRNA-seq Analysis Pipeline

Carceles-Cordon M et al., Neuron, 2026

==============================================================================

This repository contains the R scripts used for the quality control, integration, and differential gene expression analysis of single-nucleus RNA sequencing (snRNA-seq) data. The analysis focuses on Microglia from Parkinson's disease (PD) and healthy control (NC) brains across multiple studies and brain regions.

==============================================================================
1. REPOSITORY STRUCTURE
==============================================================================

The analysis assumes the following directory structure:

.
├── README.txt
├── renv.lock                          (Lockfile for exact package versions)
├── renv/                            	(Project library)
├── scripts/
│   ├── 01_qc_and_filtering.R          (QC for all raw datasets)
│   ├── 02_integration_multiregion.R   (Microglia subset, integration, and DE)
│   └── 03_integration_midbrain.R      (Midbrain-only integration and DE)
├── data/
│   ├── raw/                         (Original .mtx, .tsv, and metadata files)
│   └── processed/               (Clean Seurat objects: *_clean.rds, *_mg.rds)
└── results/
    ├── figures/                       (UMAPs, Volcano Plots, QC plots)
    └── tables/                        (Differential Expression CSV files)
==============================================================================
2. PUBLICLY AVAILABLE DATASETS USED
==============================================================================

The following public snRNA-seq datasets are analyzed in this pipeline:

1. Kamath et al. (2022)
   - Accession: GSE178265
   - Region: Midbrain 

2. Smajić et al. (2022) 
   - Accession: GSE157783
   - Region: Midbrain

3. Martirosyan et al. (2024)
   - Accession: GSE243639
   - Region: Midbrain

4. Zhang et al. (2024)
   - Accession: GSE202210
   - Region: Prefrontal cortex

==============================================================================
3. DEPENDENCIES
==============================================================================

The pipeline requires the following R packages:

- Seurat (v4.0 or higher)     : Core single-cell analysis.
- harmony (v0.1 or higher)    : Batch correction integration.
- tidyverse (v1.3 or higher)  : Data manipulation (dplyr, ggplot2).
- here (v1.0 or higher)       : Robust file path handling.
- ggrepel (v0.9 or higher)    : Non-overlapping text labels for plots.
- biomaRt (Optional)          : For gene ID conversion in Script 01.

This project uses the 'renv' package to manage R dependencies. The exact versions of all packages used in this analysis are recorded in the 'renv.lock' file, which can be accessed by using the following command in your environment:

> renv::restore() 

==============================================================================
4. ANALYSIS PIPELINE EXECUTION
==============================================================================

Scripts must be run sequentially.

------------------------------------------------------------------------------
SCRIPT 01: Quality Control and Filtering (01_qc_and_filtering.R)
------------------------------------------------------------------------------
GOAL: 
Load raw data for all four datasets (Smajic, Martirosyan, Kamath, Zhang), 
apply strict QC filtering, calculate mitochondrial content, and standardize 
metadata (Diagnosis, Donor ID, Dataset Source).

KEY STEPS:
- Load raw matrices and create initial Seurat objects.
- Convert gene identifiers (e.g., Ensembl to Symbol) where necessary.
- Calculate QC metrics (nFeature_RNA, percent.mt).
- Filter cells: nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15.

OUTPUTS:
- data/processed/*_clean.rds (One clean object per dataset)
- results/figures/qc_*_prefilter.png (QC violin plots)

------------------------------------------------------------------------------
SCRIPT 02: Multi-Region Microglia Integration (02_integration_multiregion.R)
------------------------------------------------------------------------------
GOAL: 
Subset Microglia from all four datasets, merge them, harmonize batch effects, 
and perform initial Differential Gene Expression (DGE) for PD vs NC.

INTEGRATION PARAMETERS:
- Scaling: Regress out "percent.mt".
- PCA: 30 Principal Components (PCs).
- Harmony: Run on Dims 1:30 (group.by.vars = "dataset_source").
- UMAP & Neighbors: Run on Dims 1:20.
- Clustering: Resolution 0.5.

KEY STEPS:
- Individual clustering of each dataset to manually identify Microglia.
- Subsetting and saving individual Microglia objects (*_mg.rds).
- Merging objects into 'combined_ALL_mg'.
- Running Harmony integration.
- DGE Analysis: PD vs NC (logfc.threshold = 0.05, min.pct = 0.05).

OUTPUTS:
- data/processed/combined_ALL_mg_harmony.rds
- results/figures/umap_ALL_mg_*.png (Dataset, Cluster, Region UMAPs)
- results/tables/DE_PD_vs_NC_Microglia.csv
- results/figures/Volcano_PD_vs_NC.png

------------------------------------------------------------------------------
SCRIPT 03: Midbrain-Only Microglia Integration (03_integration_midbrain.R)
------------------------------------------------------------------------------
GOAL: 
Focus exclusively on the three midbrain datasets (Smajic, Martirosyan, Kamath), 
re-harmonizing and analyzing them for region-specific validation.

PARAMETERS:
- Identical to Script 02 (30 PCs, Harmony Dims 1:30, UMAP Dims 1:20, Res 0.5).

KEY STEPS:
- Load pre-saved midbrain Microglia objects (smajic_mg, martirosyan_mg, kamath_mg).
- Merge into 'midbrain_mg'.
- Run Harmony integration workflow.
- Midbrain-specific PD vs NC DGE analysis.

OUTPUTS:
- data/processed/midbrain_mg_harmony.rds
- results/figures/Midbrain_UMAP_*.png
- results/tables/DE_Midbrain_PD_vs_NC.csv
- results/figures/Volcano_Midbrain_PD_vs_NC.png