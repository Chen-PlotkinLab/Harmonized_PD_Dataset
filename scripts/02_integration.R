# ==========================================================================
# SCRIPT 04: MICROGLIA MIDBRAIN PD DATASET HARMONIZATION
# ==========================================================================
# Goal: Load clean datasets, merge them, correct batch effects (Harmony),
#       and generate the initial UMAP visualization.

# 1. Setup and Libraries
# --------------------------------------------------------------------------
library(Seurat)
library(harmony)   # The integration tool
library(tidyverse) # For data manipulation
library(here)      # For safe file paths

options(future.globals.maxSize = 32000 * 1024^2)
message("--- Starting Script 02: Integration ---")

# 2. Load Processed Data
# --------------------------------------------------------------------------
# We load the clean objects created in Script 01.

message("Loading Smajic dataset...")
smajic <- readRDS(here("data", "processed", "smajic_clean.rds"))

message("Loading Martirosyan dataset...")
martirosyan <- readRDS(here("data", "processed", "martirosyan_clean.rds"))

message("Loading Kamath dataset...")
kamath <- readRDS(here("data", "processed", "kamath_clean.rds"))

# 3. Dataset Verification
# --------------------------------------------------------------------------
message("Verifying dimensions...")
print(paste("Smajic cells:", ncol(smajic)))
print(paste("Martirosyan cells:", ncol(martirosyan)))
print(paste("Kamath cells:", ncol(kamath)))

# Check that columns essential for integration exist in all objects
# We need: 'ID' (Sample/Donor), 'Diagnosis' (Condition), 'dataset_source' (Batch Key)
check_cols <- c("ID", "Diagnosis", "dataset_source")

for (obj_name in c("smajic", "martirosyan", "kamath")) {
  obj <- get(obj_name)
  missing <- setdiff(check_cols, colnames(obj@meta.data))
  if (length(missing) > 0) {
    stop(paste("Error:", obj_name, "is missing metadata columns:", paste(missing, collapse=", ")))
  }
}

rm(obj)
gc()

message("All datasets loaded and verified. Ready to merge.")

# 4. Merge Datasets
# --------------------------------------------------------------------------
message("Merging all three Seurat objects into 'combined' object...")
# Merges the three objects. 'add.cell.ids' prefixes the barcodes 
# (e.g., cell1 becomes Smajic_cell1) to keep track of the source.
combined <- merge(smajic, 
                  y = c(martirosyan, kamath), 
                  add.cell.ids = c("Smajic", "Martirosyan", "Kamath"), 
                  project = "PD_Midbrain_Integrated")

rm(smajic, martirosyan, kamath)
gc()

message(paste("Combined object created with", ncol(combined), "cells."))

# 5. Normalize and Find Highly Variable Features (HVGs)
# --------------------------------------------------------------------------
# Normalize data
message("Normalizing data...")
combined <- NormalizeData(combined)

# Find highly variable genes
message("Finding 2000 Highly Variable Genes (HVGs)...")
combined <- FindVariableFeatures(combined, 
                                 selection.method = "vst", 
                                 nfeatures = 2000)

# Scale Data
message("Scaling data...")
combined <- ScaleData(combined, 
                      features = VariableFeatures(combined), 
                      verbose = FALSE)

# Save highly variable gene list
hvg_list <- VariableFeatures(combined)
write_tsv(data.frame(gene = hvg_list), here("data", "processed", "integrated_HVGs.tsv"))

message("Data normalized, scaled, and HVGs identified. Ready for Harmony.")

# 6. Run Principal Component Analysis (PCA)
# --------------------------------------------------------------------------
message("Running PCA on 2000 variable features...")
combined <- RunPCA(combined, 
                   features = VariableFeatures(combined), 
                   verbose = FALSE, 
                   npcs = 50) # Using 50 Principal Components (PCs)

# 7. Run Harmony for Batch Correction
# --------------------------------------------------------------------------
message("Running Harmony integration...")
combined <- RunHarmony(combined, 
                       group.by.vars = "dataset_source", 
                       assay.use = "RNA",
                       max.iter.harmony = 20,
                       plot_convergence = FALSE)

message("Harmony integration complete. Batch effects corrected.")

# 8. Clustering and UMAP
# --------------------------------------------------------------------------

# Builds a "Shared Nearest Neighbor" (SNN) graph based on the Harmony dimensions.
message("Finding nearest neighbors using Harmony reduction...")
combined <- FindNeighbors(combined, 
                          reduction = "harmony", 
                          dims = 1:50)

message("Finding clusters at resolution 0.8...")
combined <- FindClusters(combined, 
                         resolution = 0.8)

message("Running UMAP visualization...")
combined <- RunUMAP(combined, 
                    reduction = "harmony", 
                    dims = 1:50)

message("UMAP and Clustering complete.")

# 9. Visualization & Final Save
# --------------------------------------------------------------------------

message("Generating Batch Check Plot...")
plot_batch <- DimPlot(combined, 
                      reduction = "umap", 
                      group.by = "dataset_source", 
                      label = FALSE,
                      raster = FALSE) + # raster=FALSE prevents pixelation
  ggtitle("Integration Check: UMAP by Dataset")


ggsave(here("results", "figures", "umap_batch_check.png"), 
       plot = plot_batch, 
       width = 10, 
       height = 8)

message("Generating Cluster Plot...")
plot_cluster <- DimPlot(combined, 
                        reduction = "umap", 
                        group.by = "seurat_clusters", 
                        label = TRUE, 
                        label.size = 4,
                        repel = TRUE,
                        raster = FALSE) + 
  ggtitle("UMAP Colored by Cluster (Res 0.8)")

ggsave(here("results", "figures", "umap_clusters_res08.png"), 
       plot = plot_cluster, 
       width = 10, 
       height = 8)

saveRDS(combined, file = here("data", "processed", "combined_harmony.rds"))
message("Integrated object saved to 'data/processed/combined_harmony.rds'.")

rm(combined, plot_batch, plot_cluster)
gc()

message("--- SCRIPT 02 COMPLETE: Integration and UMAP generated! ---")
