# ==========================================================================
# SCRIPT 03: MIDBRAIN MICROGLIA INTEGRATION (SMAJIC, MARTIROSYAN, KAMATH)
# ==========================================================================
# Goal: Load pre-saved Microglia objects (excluding Zhang), harmonize them,
#       generate UMAPs, and perform PD vs NC Differential Expression.

# 1. Setup and Libraries
# --------------------------------------------------------------------------
library(Seurat)
library(harmony)   
library(tidyverse) 
library(here)  
library(ggrepel)

options(future.globals.maxSize = 64000 * 1024^2) 

message("--- Starting Script 03: Midbrain Microglia Integration ---")

# 2. Load Pre-Subsetted Microglia Objects
# --------------------------------------------------------------------------
# We load the specific microglia objects saved in Script 02.

message("Loading Smajic Microglia...")
smajic_mg <- readRDS(here("data", "processed", "smajic_mg.rds"))

message("Loading Martirosyan Microglia...")
martirosyan_mg <- readRDS(here("data", "processed", "martirosyan_mg.rds"))

message("Loading Kamath Microglia...")
kamath_mg <- readRDS(here("data", "processed", "kamath_mg.rds"))

# 3. Integration (Merge & Process)
# --------------------------------------------------------------------------
message("Merging 3 Midbrain datasets into 'midbrain_mg'...")

midbrain_mg <- merge(smajic_mg, 
                     y = c(martirosyan_mg, kamath_mg), 
                     add.cell.ids = c("Smajic", "Martirosyan", "Kamath"), 
                     project = "PD_Midbrain_Microglia")

# Cleanup memory
rm(smajic_mg, martirosyan_mg, kamath_mg)
gc()

message(paste("Midbrain object created with", ncol(midbrain_mg), "cells."))

# Ensure 'dataset_source' exists
if (!"dataset_source" %in% colnames(midbrain_mg@meta.data)) {
  midbrain_mg$dataset_source <- sapply(strsplit(Cells(midbrain_mg), "_"), `[`, 1)
}

# 1. Normalize and Find HVGs
message("Normalizing data...")
midbrain_mg <- NormalizeData(midbrain_mg, verbose = FALSE)

message("Finding 2000 Highly Variable Genes...")
midbrain_mg <- FindVariableFeatures(midbrain_mg, 
                                    selection.method = "vst", 
                                    nfeatures = 2000,
                                    verbose = FALSE)

# 2. Scaling and PCA 
message("Scaling data with percent.mt regression...")
midbrain_mg <- ScaleData(midbrain_mg, 
                         features = VariableFeatures(midbrain_mg), 
                         vars.to.regress = "percent.mt", 
                         verbose = FALSE)

message("Running PCA (30 PCs)...")
midbrain_mg <- RunPCA(midbrain_mg, 
                      features = VariableFeatures(midbrain_mg), 
                      verbose = FALSE, 
                      npcs = 30)

# 3. Run Harmony (Batch Correction)
message("Running Harmony on 'dataset_source' (Dims 1:30)...")
midbrain_mg <- RunHarmony(midbrain_mg, 
                          group.by.vars = "dataset_source", 
                          assay.use = "RNA",
                          dims.use = 1:30, 
                          max.iter.harmony = 20,
                          plot_convergence = FALSE)

# 4. Clustering and UMAP 
message("Finding Neighbors (Dims 1:20)...")
midbrain_mg <- FindNeighbors(midbrain_mg, 
                             reduction = "harmony", 
                             dims = 1:20)

message("Finding Clusters (Resolution 0.5)...")
midbrain_mg <- FindClusters(midbrain_mg, 
                            resolution = 0.5, 
                            verbose = FALSE) 

message("Running UMAP (Dims 1:20)...")
midbrain_mg <- RunUMAP(midbrain_mg, 
                       reduction = "harmony", 
                       dims = 1:20, 
                       verbose = FALSE)

# 4. Visualization
# --------------------------------------------------------------------------
message("Generating UMAP Plots...")

# Plot 1: Colored by Dataset
p1 <- DimPlot(midbrain_mg, reduction = "umap", group.by = "dataset_source", raster = FALSE) + 
  ggtitle("Midbrain Integration: By Dataset")
ggsave(here("results", "figures", "Midbrain_UMAP_Dataset.png"), plot = p1, width = 10, height = 8)

# Plot 2: Colored by Clusters
p2 <- DimPlot(midbrain_mg, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE) + 
  ggtitle("Midbrain Clustering (Res 0.5)")
ggsave(here("results", "figures", "Midbrain_UMAP_Clusters.png"), plot = p2, width = 10, height = 8)

# Save the integrated object
message("Saving integrated midbrain object...")
saveRDS(midbrain_mg, file = here("data", "processed", "midbrain_mg_harmony.rds"))


# 5. Differential Expression Analysis (PD vs NC)
# --------------------------------------------------------------------------
message("--- Starting Differential Expression Analysis (PD vs NC) ---")

# 1. Prepare Metadata
# Standardize diagnosis labels to "PD" and "NC"
midbrain_mg$diagnosis_group <- case_when(
  grepl("PD|Parkinson|Disease|LBD", midbrain_mg$Diagnosis, ignore.case = TRUE) ~ "PD",
  grepl("Control|NC|Healthy", midbrain_mg$Diagnosis, ignore.case = TRUE) ~ "NC",
  TRUE ~ NA_character_
)

# Filter for PD/NC only
if (any(is.na(midbrain_mg$diagnosis_group))) {
  message("Filtering out undefined diagnoses...")
  midbrain_mg <- subset(midbrain_mg, subset = diagnosis_group %in% c("PD", "NC"))
}

# 2. Run FindMarkers
message("Setting identity to 'diagnosis_group'...")
Idents(midbrain_mg) <- "diagnosis_group"
DefaultAssay(midbrain_mg) <- "RNA"

# Check if JoinLayers is needed (Seurat v5 specific)
if (packageVersion("Seurat") >= "5.0.0") {
  message("Seurat v5 detected: Joining layers before DE...")
  midbrain_mg <- JoinLayers(midbrain_mg)
}

message("Running FindMarkers (PD vs NC)...")
de_pd_vs_nc <- FindMarkers(
  midbrain_mg,
  ident.1 = "PD",
  ident.2 = "NC",
  logfc.threshold = 0.05,
  min.pct = 0.05
)

# 3. Prepare for Volcano Plot
de_pd_vs_nc$gene <- rownames(de_pd_vs_nc)
de_pd_vs_nc$log10_pval <- -log10(de_pd_vs_nc$p_val_adj)

# Define genes to highlight (GPNMB + others from Harmony_All)
highlight_genes <- c(
  "GPNMB",  
  "CX3CR1", "P2RY12")

de_pd_vs_nc$highlight <- ifelse(de_pd_vs_nc$gene %in% highlight_genes, "yes", "no")

# Save table
write.csv(de_pd_vs_nc, here("results", "tables", "DE_Midbrain_PD_vs_NC.csv"))

# 4. Generate Volcano Plot
message("Generating Volcano Plot...")

p_volcano <- ggplot(de_pd_vs_nc, aes(x = avg_log2FC, y = log10_pval)) +
  geom_point(aes(color = highlight), alpha = 0.7) +
  scale_color_manual(values = c("no" = "gray", "yes" = "red")) +
  geom_text_repel(
    data = subset(de_pd_vs_nc, highlight == "yes"),
    aes(label = gene),
    size = 4, 
    box.padding = 0.3, 
    max.overlaps = 20
  ) +
  theme_minimal() +
  labs(
    title = "PD vs NC: Midbrain Microglia",
    x = "Average log2 Fold Change",
    y = "-log10 Adjusted p-value"
  )

ggsave(here("results", "figures", "Volcano_Midbrain_PD_vs_NC.png"), 
       plot = p_volcano, 
       width = 8, height = 6, dpi = 300)

message("--- SCRIPT 03 COMPLETE: Midbrain Analysis Finished ---")
