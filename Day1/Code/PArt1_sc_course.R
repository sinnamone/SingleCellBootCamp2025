# Load required libraries
library(Seurat)
library(dplyr)

# ----------------------------
# STEP 1: Load PBMC 3k dataset
# ----------------------------
pbmc.data <- Read10X(data.dir = "C:/Users/remoteadmin/Downloads/wetransfer_filtered_gene_bc_matrices_2025-07-14_1142/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "PBMC3K", min.cells = 3, min.features = 200)

# ----------------------------
# STEP 2: Basic object summary
# ----------------------------

# Print general information about the Seurat object
pbmc

# ----------------------------
# STEP 3: Explore cell and gene names
# ----------------------------

# Print first 10 cell barcodes (column names = cells)
head(colnames(pbmc), 10)

# Print first 10 gene names (row names = genes)
head(rownames(pbmc), 10)

# ----------------------------
# STEP 4: Explore metadata
# ----------------------------

# Show the first rows of the metadata (one row per cell)
head(pbmc@meta.data)

# Add a new metadata column (e.g., random group labels)
set.seed(42)  # for reproducibility
pbmc$group <- sample(c("GroupA", "GroupB"), ncol(pbmc), replace = TRUE)

# View the new column
head(pbmc$group)

# ----------------------------
# STEP 5: Explore assays and expression data
# ----------------------------

# List all assays (typically just "RNA" at this stage)
Assays(pbmc)

# View raw counts (gene expression before normalization)
GetAssayData(pbmc, slot = "counts")[1:5, 1:5]        # ✅ Seurat v5
# pbmc[["RNA"]]@counts[1:5, 1:5]                      # ❌ Seurat v4 style – won't work in v5

# (Optional: run normalization)
pbmc <- NormalizeData(pbmc)

# View normalized data
GetAssayData(pbmc, slot = "data")[1:5, 1:5]           # ✅ Seurat v5
# pbmc[["RNA"]]@data[1:5, 1:5]                        # ❌ Seurat v4 style – won't work in v5

# ----------------------------
# STEP 6: Feature and reduction information
# ----------------------------

# Compute variable genes
pbmc <- FindVariableFeatures(pbmc)

# Print top 10 variable features (genes)
head(VariableFeatures(pbmc), 10)

# Run PCA
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)

# List dimensional reductions (e.g., "pca")
Reductions(pbmc)

# ----------------------------
# STEP 7: Idents and clusters
# ----------------------------

# Print current identities (clusters, if defined)
Idents(pbmc)

# Check slot names available in Seurat v5 assay
slotNames(pbmc[["RNA"]])

# Check available assay data slots (valid for Seurat v5):
GetAssay(pbmc, assay = "RNA")  # prints a summary of the RNA assay

GetAssayData(pbmc, assay = "RNA", slot = "counts")      # raw counts
GetAssayData(pbmc, assay = "RNA", slot = "data")        # normalized data
GetAssayData(pbmc, assay = "RNA", slot = "scale.data")  # scaled data

