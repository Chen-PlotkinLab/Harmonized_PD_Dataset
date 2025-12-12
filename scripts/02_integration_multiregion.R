# ==========================================================================
# SCRIPT 02: MULTI-REGION DATASET INTEGRATION
# ==========================================================================
# Goal: Load clean datasets, SUBSET TO MICROGLIA, merge them, correct batch 
#       effects (Harmony), and generate UMAP visualization / downstream analysis.
#       Output object: 'combined_ALL_mg'

# 1. Setup and Libraries
# --------------------------------------------------------------------------
library(Seurat)
library(harmony)   
library(tidyverse) 
library(here)  
library(ggrepel)
options(future.globals.maxSize = 64000 * 1024^2) 

message("--- Starting Script 02B: Multi-Region Microglia Integration ---")

# 2. Load Processed Data (All 4 Datasets)
# --------------------------------------------------------------------------
# Load the full clean objects 
smajic <- readRDS(here("data", "processed", "smajic_clean.rds"))
martirosyan <- readRDS(here("data", "processed", "martirosyan_clean.rds"))
kamath <- readRDS(here("data", "processed", "kamath_clean.rds"))
zhang <- readRDS(here("data", "processed", "zhang_clean.rds"))

# 3. Individual Clustering & Manual Annotation (Pre-Integration)
# --------------------------------------------------------------------------
# Goal: Cluster each dataset individually, plot markers, and manually select Microglia.

# A. Define Broad Markers 
markers <- list(
  "Microglia"        = c("CSF1R", "P2RY12", "TMEM119", "CX3CR1", "CD74"),
  "Neuron"           = c("SYT1", "SNAP25", "GRIN1"),
  "Astrocyte"        = c("AQP4", "GFAP", "SLC1A3"),
  "Oligodendrocyte"  = c("MBP", "MOG", "PLP1", "MOBP"),
  "Endothelial"      = c("CLDN5", "VWF", "PECAM1", "FLT1")
)

# B. Helper Function: Process, Cluster, and Plot
generate_diagnostic_plots <- function(obj, name, n_pcs = 15, res = 0.5) {
  
  message(paste0(">>> Processing ", name, "..."))
  
  # Standard Workflow
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:n_pcs, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res, verbose = FALSE) 
  obj <- RunUMAP(obj, dims = 1:n_pcs, verbose = FALSE)
  
  # 1. DotPlot (To identify cell types)
  p_dot <- DotPlot(obj, features = markers, group.by = "seurat_clusters") + 
    RotatedAxis() +
    ggtitle(paste(name, "- Marker Expression"))
  
  # 2. UMAP (To see cluster structure)
  p_umap <- DimPlot(obj, label = TRUE, label.size = 4) + 
    ggtitle(paste(name, "- Clusters"))
  
  # Save plots
  ggsave(here("results", "figures", paste0("QC_Annotation_", name, ".png")), 
         plot = p_dot + p_umap, 
         width = 14, height = 7)
  
  message(paste0("    Plots saved to: results/figures/QC_Annotation_", name, ".png"))
  return(obj)
}

# C. Run Processing on All 4 Datasets
# Note: This overwrites the objects with their clustered versions
smajic      <- generate_diagnostic_plots(smajic, "Smajic")
martirosyan <- generate_diagnostic_plots(martirosyan, "Martirosyan")
kamath      <- generate_diagnostic_plots(kamath, "Kamath")
zhang       <- generate_diagnostic_plots(zhang, "Zhang")

message("\n========================================================")
message("🛑 STOP! CHECK YOUR PLOTS NOW.")
message("1. Go to 'results/figures/' folder.")
message("2. Open the 4 'QC_Annotation_*.png' files.")
message("3. Identify which cluster numbers correspond to 'Microglia'.")
message("4. Fill in the 'ids_to_keep' variables below.")
message("========================================================\n")


# 4. Manual Selection 
# ----------------------------------------------------------------------
# Enter the cluster numbers inside c(). Example: c(0, 1, 4)

ids_smajic      <- c(3, 10)
ids_martirosyan <- c(17, 7, 4)  
ids_kamath      <- c(2, 7, 21)
ids_zhang       <- c(16, 5, 3, 20)

message("Subsetting datasets to selected Microglia clusters...")

smajic_mg      <- subset(smajic, idents = ids_smajic)
martirosyan_mg <- subset(martirosyan, idents = ids_martirosyan)
kamath_mg      <- subset(kamath, idents = ids_kamath)
zhang_mg       <- subset(zhang, idents = ids_zhang)

# Tag them for the integration step
smajic_mg$celltype_group      <- "Microglia"
martirosyan_mg$celltype_group <- "Microglia"
kamath_mg$celltype_group      <- "Microglia"
zhang_mg$celltype_group       <- "Microglia"

# Cleanup memory
rm(smajic, martirosyan, kamath, zhang)
gc()

message("Saving individual Microglia objects before integration...")

saveRDS(smajic_mg,      file = here("data", "processed", "smajic_mg.rds"))
saveRDS(martirosyan_mg, file = here("data", "processed", "martirosyan_mg.rds"))
saveRDS(kamath_mg,      file = here("data", "processed", "kamath_mg.rds"))
saveRDS(zhang_mg,       file = here("data", "processed", "zhang_mg.rds"))

message("Individual microglia objects saved.")

# 5. Integration
# --------------------------------------------------------------------------
message("Merging all 4 Microglia datasets into 'combined_ALL_mg'...")

combined_all_mg <- merge(smajic_mg, 
                         y = c(martirosyan_mg, kamath_mg, zhang_mg), 
                         add.cell.ids = c("Smajic", "Martirosyan", "Kamath", "Zhang"), 
                         project = "PD_MultiRegion_Microglia")

# Cleanup memory
rm(smajic_mg, martirosyan_mg, kamath_mg, zhang_mg)
gc()

message(paste("combined_ALL_mg object created with", ncol(combined_all_mg), "cells."))

# Ensure 'dataset_source' exists for Harmony (extracting from prefix if needed)
if (!"dataset_source" %in% colnames(combined_all_mg@meta.data)) {
  combined_all_mg$dataset_source <- sapply(strsplit(Cells(combined_all_mg), "_"), `[`, 1)
}

# 1. Normalize and Find Highly Variable Features (HVGs)
message("Normalizing data...")
combined_all_mg <- NormalizeData(combined_all_mg, verbose = FALSE)

message("Finding 2000 Highly Variable Genes (HVGs)...")
combined_all_mg <- FindVariableFeatures(combined_all_mg, 
                                        selection.method = "vst", 
                                        nfeatures = 2000,
                                        verbose = FALSE)

# 2. Scaling and PCA 
message("Scaling data with percent.mt regression...")
combined_all_mg <- ScaleData(combined_all_mg, 
                             features = VariableFeatures(combined_all_mg), 
                             vars.to.regress = "percent.mt", 
                             verbose = FALSE)

message("Running PCA (30 PCs)...")
combined_all_mg <- RunPCA(combined_all_mg, 
                          features = VariableFeatures(combined_all_mg), 
                          verbose = FALSE, 
                          npcs = 30)

# 3. Run Harmony (Batch Correction)
message("Running Harmony integration on 'dataset_source' (Dims 1:30)...")
combined_all_mg <- RunHarmony(combined_all_mg, 
                              group.by.vars = "dataset_source", 
                              assay.use = "RNA",
                              dims.use = 1:30, 
                              max.iter.harmony = 20,
                              plot_convergence = FALSE)

message("Harmony complete.")

# 4. Clustering and UMAP (Matching Harmony_All.R)
message("Finding Neighbors (Dims 1:20)...")
combined_all_mg <- FindNeighbors(combined_all_mg, 
                                 reduction = "harmony", 
                                 dims = 1:20)

message("Finding Clusters (Resolution 0.5)...")
combined_all_mg <- FindClusters(combined_all_mg, 
                                resolution = 0.5, 
                                verbose = FALSE) 

message("Running UMAP (Dims 1:20)...")
combined_all_mg <- RunUMAP(combined_all_mg, 
                           reduction = "harmony", 
                           dims = 1:20, 
                           verbose = FALSE)

# 5. Visualization & Save

# Plot 1: Colored by Dataset
p1 <- DimPlot(combined_all_mg, reduction = "umap", group.by = "dataset_source", raster = FALSE) + 
  ggtitle("Microglia Integration Check: By Dataset")
ggsave(here("results", "figures", "umap_ALL_mg_dataset.png"), plot = p1, width = 10, height = 8)

# Plot 2: Colored by Clusters
p2 <- DimPlot(combined_all_mg, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE) + 
  ggtitle("Multi-Region Microglia Clustering (Res 0.5)")
ggsave(here("results", "figures", "umap_ALL_mg_clusters.png"), plot = p2, width = 10, height = 8)

# Save Object
message("Saving combined_ALL_mg object...")
saveRDS(combined_all_mg, file = here("data", "processed", "combined_ALL_mg_harmony.rds"))

message("--- SCRIPT 02 COMPLETE ---")


# 6. Differential Expression Analysis
# --------------------------------------------------------------------------

message("--- Starting Differential Expression Analysis ---")

# 1. Prepare metadata

combined_all_mg$diagnosis_group <- case_when(
  grepl("PD|Parkinson|Disease|LBD", combined_all_mg$Diagnosis, ignore.case = TRUE) ~ "PD",
  grepl("Control|NC|Healthy", combined_all_mg$Diagnosis, ignore.case = TRUE) ~ "NC",
  TRUE ~ NA_character_
)

# Remove any cells that couldn't be classified (optional safety step)
if (any(is.na(combined_all_mg$diagnosis_group))) {
  message("Warning: Some cells has undefined diagnosis. Subsetting to PD/NC only.")
  combined_all_mg <- subset(combined_all_mg, subset = diagnosis_group %in% c("PD", "NC"))
}

# 2. Run FindMarkers

message("Setting identity to 'diagnosis_group'...")
Idents(combined_all_mg) <- "diagnosis_group"

# Ensure we use the RNA assay and join layers
DefaultAssay(combined_all_mg) <- "RNA"

# Check if JoinLayers is needed (Seurat v5 specific)
if (packageVersion("Seurat") >= "5.0.0") {
  message("Seurat v5 detected: Joining layers before DE...")
  combined_all_mg <- JoinLayers(combined_all_mg)
}

message("Running FindMarkers (PD vs NC)...")
# Note: Low thresholds (0.05) as per your script
de_pd_vs_nc <- FindMarkers(
  combined_all_mg,
  ident.1 = "PD",
  ident.2 = "NC",
  logfc.threshold = 0.05,
  min.pct = 0.05
)

# 3. DGE Table

de_pd_vs_nc$gene <- rownames(de_pd_vs_nc)
de_pd_vs_nc$log10_pval <- -log10(de_pd_vs_nc$p_val_adj)
de_pd_vs_nc$highlight <- ifelse(de_pd_vs_nc$gene %in% highlight_genes, "yes", "no")

# Save the DE table
write.csv(de_pd_vs_nc, here("results", "tables", "DE_PD_vs_NC_Microglia.csv"))

# 4. Volcano Plot

message("Generating Volcano Plot...")

highlight_genes <- c(
  "GPNMB", "P2RY12", "CX3CR1"
)

de_pd_vs_nc$highlight <- ifelse(de_pd_vs_nc$gene %in% highlight_genes, "yes", "no")

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
    title = "PD vs NC Differential Expression (Microglia)",
    x = "Average log2 Fold Change",
    y = "-log10 Adjusted p-value"
  )

ggsave(here("results", "figures", "Volcano_PD_vs_NC.png"), 
       plot = p_volcano, 
       width = 8, height = 6, dpi = 300)

message("Differential Expression Analysis Complete. Results saved.")






