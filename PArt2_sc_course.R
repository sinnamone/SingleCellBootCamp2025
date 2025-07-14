# Load required libraries
library(Seurat)
library(dplyr)

# ----------------------------
# STEP 1: Load PBMC 3k dataset
# ----------------------------
pbmc.data <- Read10X(data.dir = "C:/Users/remoteadmin/Downloads/wetransfer_filtered_gene_bc_matrices_2025-07-14_1142/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "PBMC3K", min.cells = 3, min.features = 200)

# ----------------------------------------
# STEP 1: Calculate Quality Control Metrics
# ----------------------------------------

# Mitochondrial genes usually start with "MT-"
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# View updated metadata (first 6 rows)
head(pbmc@meta.data)

# ----------------------------------------
# STEP 2: Visualize QC Metrics
# ----------------------------------------

# Violin plots for key metrics: nFeature_RNA, nCount_RNA, percent.mt
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0.1)

# Scatter plots to explore potential filtering thresholds
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# ----------------------------------------
# STEP 3: Filter Low-Quality Cells
# ----------------------------------------

# Filter cells based on:
# - Min features: >200
# - Max features: <2500
# - Max mitochondrial content: <5%
pbmc_filtered <- subset(pbmc, subset = nFeature_RNA > 200 &
                          nFeature_RNA < 2500 &
                          percent.mt < 5)

# Print a summary of the filtered object
pbmc_filtered

# ----------------------------------------
# STEP 4: (Optional) Filter Lowly Expressed Genes
# ----------------------------------------

# Remove genes expressed in fewer than X cells (e.g., fewer than 3 cells)
# Note: This is not always necessary after Seurat’s default handling,
# but useful to understand filtering principles

counts_matrix <- GetAssayData(pbmc_filtered, slot = "counts")
genes_to_keep <- rowSums(counts_matrix > 0) >= 10
pbmc_filtered <- pbmc_filtered[genes_to_keep, ]

# Confirm number of genes and cells after filtering
dim(pbmc_filtered)

# ----------------------------------------
# STEP 5: Save filtered object (optional)
# ----------------------------------------

# Save for later use
saveRDS(pbmc_filtered, file = "C:/Users/remoteadmin/Downloads/pbmc_filtered_seurat.rds")
