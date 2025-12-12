## --- 01_qc_and_filtering.R ---
## Goal: QC and standardize Smajic, Martirosyan, and Kamath datasets individually.
## Note: All file paths use here::here() and assume files are correctly placed in a 'data/raw' folder.

# 0. SETUP
# --------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(here)
library(biomaRt)
library(Matrix)
library(readr)
library(readxl)

# Ensure output folders exist
if (!dir.exists(here("data", "processed"))) {
  dir.create(here("data", "processed"))
}
if (!dir.exists(here("results", "figures"))) {
  dir.create(here("results", "figures"))
}

# ==========================================================================
# DATASET 1: SMAJIC (GSE157783) 
# ==========================================================================
message("--- Starting Smajic (GSE157783) Dataset Processing ---")

# A. Load Data, Convert Genes, and Create Seurat Object
# --------------------------------------------------------------------------
smajic_data_dir <- here("data", "raw", "GSE157783_Smajic")

# 1. Load the raw files
message("Loading raw TSV files...")
barcode_meta <- readr::read_delim(here(smajic_data_dir, "IPDCO_hg_midbrain_cell.tsv"), 
                                  delim = "\t", show_col_types = FALSE)
genes <- readr::read_delim(here(smajic_data_dir, "IPDCO_hg_midbrain_genes.tsv"), 
                           delim = "\t", col_names = FALSE, show_col_types = FALSE)
umi_mat_df <- readr::read_delim(here(smajic_data_dir, "IPDCO_hg_midbrain_UMI.tsv"), 
                                delim = "\t", show_col_types = FALSE)

# 2. Fix Matrix Orientation
gene_ids <- make.unique(genes$X1)[-1] 

# Transpose double-shuffle to get (Genes x Cells)
umi_mat <- as.matrix(t(umi_mat_df)) 
colnames(umi_mat) <- gene_ids
rownames(umi_mat) <- colnames(umi_mat_df)
mat <- t(umi_mat) 
mat <- as(mat, "dgCMatrix")

# 3. Gene Symbol Conversion
all_genes <- rownames(mat)
ensg_indices <- grep("^ENSG", all_genes)
ensg_ids_to_convert <- all_genes[ensg_indices]
ensg_base <- gsub("\\..*$", "", ensg_ids_to_convert) 

message(paste("Identified", length(ensg_indices), "Ensembl IDs to convert..."))

if (length(ensg_indices) > 0) {
  # Connect to Ensembl
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  # Query ONLY the Ensembl IDs
  gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = ensg_base,
                    mart = mart)
  
  # Create lookup table
  symbol_lookup <- setNames(gene_map$hgnc_symbol, gene_map$ensembl_gene_id)
  
  # Apply conversion to the matrix row names
  new_gene_names <- all_genes
  
  for (i in seq_along(ensg_indices)) {
    idx <- ensg_indices[i]
    id_base <- ensg_base[i]
    
    if (id_base %in% names(symbol_lookup)) {
      sym <- symbol_lookup[[id_base]]
      if (!is.na(sym) && sym != "") {
        new_gene_names[idx] <- sym
      }
    }
  }
  
  # Handle duplicates (e.g. multiple ENSG IDs mapping to same symbol)
  rownames(mat) <- make.unique(new_gene_names)
}

# 4. Create Seurat Object (Using original filters)
barcode_meta_df <- as.data.frame(barcode_meta)
rownames(barcode_meta_df) <- barcode_meta_df$barcode
barcode_meta_final <- barcode_meta_df[colnames(mat), ]

smajic <- CreateSeuratObject(
  counts = mat, 
  project = "Smajic", 
  min.cells = 3,        
  min.features = 200,   
  meta.data = barcode_meta_final
)

message("Seurat Object created successfully.")

# Cleanup
rm(barcode_meta_df, barcode_meta_final, umi_mat_df, umi_mat, genes, gene_ids, mat, mart, gene_map)
gc()





# B. Standardized Metadata Merge 
# --------------------------------------------------------------------------
sample_info_path <- here("data", "raw", "GSE157783_Smajic", "Sample_Info.csv")
metadata_df <- read_csv(sample_info_path, show_col_types = FALSE) %>%
  dplyr::select(ID, Diagnosis) %>% 
  distinct()

# Join using 'patient' (from original file) which will be substituted for 'ID' in clean metadata
meta_to_join <- smajic@meta.data %>%
  dplyr::select(patient) %>%           
  tibble::rownames_to_column("barcode") %>%
  left_join(metadata_df, by = c("patient" = "ID")) %>% 
  tibble::column_to_rownames("barcode")
smajic$Diagnosis <- meta_to_join$Diagnosis
smajic$ID <- smajic$patient            # Rename 'patient' to standard 'ID' 
smajic$dataset_source <- "Smajic"
smajic$Diagnosis <- factor(smajic$Diagnosis, levels = c("NC", "PD"))
smajic$patient <- NULL

message("Smajic data loaded and metadata merged successfully.")

# C. QC & Filtering 
# --------------------------------------------------------------------------
# 1. Calculate QC metrics

smajic[["percent.mt"]] <- PercentageFeatureSet(smajic, pattern = "^MT-")

plot_qc <- VlnPlot(smajic, 
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                   ncol = 3, 
                   group.by = "Diagnosis", 
                   pt.size = 0.01)
ggsave(here("results", "figures", "qc_smajic_prefilter.png"), plot = plot_qc, width = 10, height = 5)

# 2. Filter
message(paste("Cells before filtering:", ncol(smajic)))

smajic <- subset(smajic, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 3000 &  # We chose to keep this parameter consistent, but can be changed
                   percent.mt < 15)       

message(paste("Cells after filtering:", ncol(smajic)))

# 3. Save
saveRDS(smajic, file = here("data", "processed", "smajic_clean.rds"))
message("Smajic dataset processed and saved to 'data/processed/smajic_clean.rds'.")
rm(smajic, plot_qc)
gc()


# ==========================================================================
# DATASET 2: MARTIROSYAN (GSE243639)
# ==========================================================================
message("--- Starting Martirosyan (GSE243639) Dataset Processing ---")

# A. Automatic Organization, Loading, and Merging into Seurat Object
# --------------------------------------------------------------------------
martirosyan_dir <- here("data", "raw", "GSE243639_Martirosyan")

# Check and Organize Files, if needed, in the root folder.
flat_files <- list.files(martirosyan_dir, pattern = "_matrix.mtx.gz$", full.names = FALSE)

if (length(flat_files) > 0) {
  message(paste("Found", length(flat_files), "unorganized files. Reorganizing into folders..."))

  samples_to_organize <- unique(sub("_matrix.mtx.gz", "", flat_files))
  
  for (sample in samples_to_organize) {

    new_folder <- file.path(martirosyan_dir, sample)
    if (!dir.exists(new_folder)) dir.create(new_folder)
    
    files_to_move <- c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
    
    for (f_type in files_to_move) {
      original_file <- paste0(sample, "_", f_type)
      source_path <- file.path(martirosyan_dir, original_file)
      dest_path   <- file.path(new_folder, f_type) 

      if (file.exists(source_path)) {
        file.rename(from = source_path, to = dest_path)
      }
    }
  }
  message("Organization complete. Proceeding to load.")
} else {
  message("Files appear already organized. Proceeding to load.")
}

# Now we load the data in 10X format
sample_folders <- list.dirs(martirosyan_dir, recursive = FALSE, full.names = TRUE)
is_sample_folder <- sapply(sample_folders, function(x) file.exists(file.path(x, "matrix.mtx.gz")))
sample_folders <- sample_folders[is_sample_folder]

if (length(sample_folders) == 0) {
  stop("No valid sample folders found in 'data/raw/GSE243639_Martirosyan'.")
}

seurat_list <- list()

for (folder in sample_folders) {
  sample_id <- basename(folder)
  
  tryCatch({
    raw_counts <- Read10X(data.dir = folder)
    
    # Create Object with initial filters (Matched to other datasets!)
    obj <- CreateSeuratObject(
      counts = raw_counts, 
      project = sample_id, 
      min.cells = 3, 
      min.features = 200
    )
    
    seurat_list[[sample_id]] <- obj
    
  }, error = function(e) {
    warning(paste("Could not load sample:", sample_id, "-", e$message))
  })
}

# Now we merge into a Seurat object
if (length(seurat_list) > 0) {
  message("Merging individual Martirosyan samples...")
  
  martirosyan <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list), 
    project = "Martirosyan"
  )
  
  # Create a standard 'ID' column (currently holds GSM IDs) to match other datasets
  martirosyan$ID <- martirosyan$orig.ident
  
  message(paste("Martirosyan Object Created. Cells:", ncol(martirosyan)))
  
  # Cleanup
  rm(seurat_list, raw_counts, obj)
  gc()
  
} else {
  stop("No Seurat objects were created. Check data path.")
}


# B: Standardized Metadata Merge (Martirosyan)
# --------------------------------------------------------------------------
clinical_path <- here("data", "raw", "GSE243639_Martirosyan", "GSE243639_Clinical_data.csv")

# 1. Load Clinical Data
# The file uses semicolons ';' and has ~5 lines of description at the top.
# We need to pull the sample ID from the original files 

martirosyan$ID_raw <- martirosyan$ID 
martirosyan$ID <- stringr::str_extract(martirosyan$ID_raw, "s_\\d+")

metadata_raw <- read_delim(
  clinical_path, 
  delim = ";", 
  skip = 5, 
  show_col_types = FALSE
)

message("Clinical file loaded. Columns found:")
print(colnames(metadata_raw)) 

metadata_clean <- metadata_raw %>%
  # Select the two columns needed: 'Sample ID' from the CSV and 'Clinical diagnosis'
  dplyr::select(
    ID = `Sample ID`,              # Map CSV 'Sample ID' (s.xxxx) to standard 'ID'
    Diagnosis = `Clinical diagnosis` # Map CSV 'Clinical diagnosis'
  ) %>%
  mutate(ID = stringr::str_replace(ID, "\\.", "_")) %>%
  mutate(
    # Standardize Diagnosis: 'Control' -> 'NC', 'Parkinson's' -> 'PD'
    Diagnosis = case_when(
      Diagnosis == "Control" ~ "NC",
      Diagnosis == "Parkinson's" ~ "PD",
      TRUE ~ "Other"
    )
  ) %>%
  distinct()

# 3. Merge with Seurat object
meta_to_join <- martirosyan@meta.data %>%
  dplyr::select(ID) %>% 
  tibble::rownames_to_column("barcode") %>%
  left_join(metadata_clean, by = "ID") %>% 
  tibble::column_to_rownames("barcode")

martirosyan$Diagnosis <- meta_to_join$Diagnosis
martirosyan$dataset_source <- "Martirosyan"
martirosyan$Diagnosis <- factor(martirosyan$Diagnosis, levels = c("NC", "PD"))

message("Martirosyan metadata merged.")
message("Check Diagnosis counts:")
print(table(martirosyan$Diagnosis, useNA = "ifany"))



# C. QC and filtering
# --------------------------------------------------------------------------
# 1. Calculate and visualize QC metrics

[["percent.mt"]] <- PercentageFeatureSet(martirosyan, pattern = "^MT-")
message("Generating QC Violin Plots (pre-filter)...")
plot_qc_martirosyan <- VlnPlot(martirosyan, 
                               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                               ncol = 3, 
                               group.by = "Diagnosis", 
                               pt.size = 0.01)
ggsave(here("results", "figures", "qc_martirosyan_prefilter.png"), 
       plot = plot_qc_martirosyan, 
       width = 10, 
       height = 5)

# 2. Apply filtering for harmonization (nFeature_RNA >200 and <3000; percent.mt < 15)

message(paste("Martirosyan cells before filtering:", ncol(martirosyan)))

martirosyan <- subset(martirosyan, subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 3000 &  
                        percent.mt < 15)       

message(paste("Martirosyan cells after filtering:", ncol(martirosyan)))

# 3. Save object 
saveRDS(martirosyan, file = here("data", "processed", "martirosyan_clean.rds"))
message("Martirosyan dataset processed and saved to 'data/processed/martirosyan_clean.rds'.")

rm(
  martirosyan,          # The heavy Seurat object (it's safe on disk now)
  plot_qc_martirosyan,  # The plot object
  metadata_raw,         # The raw clinical CSV data
  metadata_clean,       # The cleaned clinical data
  meta_to_join,         # The temporary joining table
  clinical_path,        # Path variables
  martirosyan_dir, 
  sample_folders
)
gc()
message("Memory cleared. Ready for the next dataset.")






# ==========================================================================
# DATASET 3: KAMATH (GSE178265)
# ==========================================================================
message("--- Starting Kamath (GSE178265) Dataset Processing ---")

# A. Load Data and Create Seurat Object
# --------------------------------------------------------------------------
kamath_dir <- here("data", "raw", "GSE178265_Kamath")

# 1. Define File Paths (Based on your list.files output)

matrix_path   <- here(kamath_dir, "GSE178265_Homo_matrix.mtx.gz")
features_path <- here(kamath_dir, "GSE178265_Homo_features.tsv.gz")
barcodes_path <- here(kamath_dir, "GSE178265_Homo_bcd.tsv.gz")

message("Loading Kamath matrix files manually...")

# 2. Read Components

mat <- Matrix::readMM(matrix_path)

# Load features (genes). Assumes standard 10x format (V1=ID, V2=Symbol)
features <- readr::read_tsv(features_path, col_names = FALSE, show_col_types = FALSE)

# Load barcodes. 
barcodes <- readr::read_tsv(barcodes_path, col_names = FALSE, show_col_types = FALSE)

# 3. Assign Row/Col Names

# Genes: Use column X2 (Symbols)
rownames(mat) <- make.unique(features$X2)

# Cells: Use column X1 (Barcodes)
colnames(mat) <- barcodes$X1

# 4. Create Seurat Object

kamath <- CreateSeuratObject(
  counts = mat,
  project = "Kamath",
  min.cells = 3,
  min.features = 200
)

message(paste("Kamath Object Created. Cells:", ncol(kamath), "Features:", nrow(kamath)))

# Cleanup large temp files
rm(mat, features, barcodes, gene_names)
gc()

# B. Load Data and Create Seurat Object
# -------------------------------------------------------------------------
# 1. Load clinical data

clinical_path <- here("data", "raw", "GSE178265_Kamath", "Sample_Info.xlsx")

metadata_raw <- readxl::read_excel(clinical_path)

metadata_clean <- metadata_raw %>%
  dplyr::select(
    ID = ID,                       # Donor ID (e.g., 4340)
    Diagnosis = Disease,           # Disease status to match with harmonized
    Sex = Sex,                     
    Age = Age                      
  ) %>%
  mutate(
    # Standardize Diagnosis
    Diagnosis = case_when(
      Diagnosis == "Control" ~ "NC",
      Diagnosis %in% c("PD", "LBD") ~ "PD", 
      TRUE ~ "Other"
    ),
    # Convert numeric ID (e.g., 4340) to character for matching
    ID = as.character(ID) 
  ) %>%
  distinct()

# 2. Merge metadata with Seurat

message("Extracting Donor ID from cell barcodes...")

# Extract the FIRST 4-digit number found in the string (e.g., 4340)
kamath$ID_extracted <- stringr::str_extract(rownames(kamath@meta.data), "\\d{4}")

message(paste("Example ID extraction:", rownames(kamath@meta.data)[1], "->", kamath$ID_extracted[1]))

meta_to_join <- kamath@meta.data %>%
  dplyr::select(ID_extracted) %>% 
  tibble::rownames_to_column("barcode") %>%
  # Join extracted ID with clinical ID
  left_join(metadata_clean, by = c("ID_extracted" = "ID")) %>% 
  tibble::column_to_rownames("barcode")

kamath$Diagnosis <- meta_to_join$Diagnosis
kamath$Sex <- meta_to_join$Sex
kamath$Age <- meta_to_join$Age
kamath$dataset_source <- "Kamath"
kamath$ID <- kamath$ID_extracted 
kamath$ID_extracted <- NULL # Cleanup

kamath$Diagnosis <- factor(kamath$Diagnosis, levels = c("NC", "PD"))

message("Kamath metadata merged.")
message("Check Diagnosis counts (Should show NC and PD):")
print(table(kamath$Diagnosis, useNA = "ifany"))

# We found 4 donors (4340, 3839, 5730, 3898) missing from the clinical metadata, which will be removed
message(paste("Total cells before metadata pruning:", ncol(kamath)))
kamath <- subset(kamath, subset = !is.na(Diagnosis))
message(paste("Total cells after removing missing donors:", ncol(kamath)))
print(table(kamath$Diagnosis, useNA = "ifany"))

# 3. Final Cleanup 
rm(metadata_raw, metadata_clean, meta_to_join, clinical_path)
gc()



# C. QC and Filtering
# -------------------------------------------------------------------------
# 1. Calculate & Visualize QC metrics 

kamath[["percent.mt"]] <- PercentageFeatureSet(kamath, pattern = "^MT-")

message("Generating QC Violin Plots (pre-filter)...")

plot_qc_kamath <- VlnPlot(kamath, 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                          ncol = 3, 
                          group.by = "Diagnosis", 
                          pt.size = 0) # Set pt.size=0 for large datasets to speed up plotting

ggsave(here("results", "figures", "qc_kamath_prefilter.png"), 
       plot = plot_qc_kamath, 
       width = 10, 
       height = 5)

# 2. Apply filtering step 

message(paste("Kamath cells before QC filtering:", ncol(kamath)))

kamath <- subset(kamath, subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 3000 &  
                   percent.mt < 15)       

message(paste("Kamath cells after QC filtering:", ncol(kamath)))

# 3. Save & Clear Memory

saveRDS(kamath, file = here("data", "processed", "kamath_clean.rds"))
message("Kamath dataset processed and saved to 'data/processed/kamath_clean.rds'.")

rm(kamath, plot_qc_kamath)
gc()

message("--- SCRIPT 01 COMPLETE: All 3 datasets processed! ---")

