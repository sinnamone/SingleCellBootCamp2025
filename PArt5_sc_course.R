# Load required libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(harmony)

# Setup options
options(timeout = 1000)
InstallData("pbmcsca")
options(future.globals.maxSize = 1e9)

# STEP 1: Load and Subset to 3 Methods
obj <- LoadData("pbmcsca")
obj <- subset(obj, nFeature_RNA > 1000)

# Select 3 available methods that are known to work
methods_to_use <- c("Smart-seq2", "10x Chromium (v3)", "CEL-Seq2")
obj_subset <- subset(obj, subset = Method %in% methods_to_use)

# Merge directly without filtering individual subsets
obj_merged <- obj_subset

# -----------------------------------------------
# STEP 2: QC and Normalization on Merged Object
# -----------------------------------------------
obj_merged[["percent.mt"]] <- PercentageFeatureSet(obj_merged, pattern = "^MT-")
obj_merged <- subset(obj_merged, subset = nFeature_RNA > 200 & percent.mt < 20)
obj_merged <- NormalizeData(obj_merged)
obj_merged <- FindVariableFeatures(obj_merged, selection.method = "vst", nfeatures = 2000)
obj_merged <- ScaleData(obj_merged)

# PCA and UMAP before integration
obj_merged <- RunPCA(obj_merged, npcs = 30)
obj_merged <- RunUMAP(obj_merged, reduction = "pca", dims = 1:20)

p_before <- DimPlot(obj_merged, reduction = "umap", group.by = "Method") + 
  ggtitle("UMAP Before Harmony Integration")

# -----------------------------------------------
# STEP 3: Harmony Integration
# -----------------------------------------------
library(harmony)
obj_merged <- RunHarmony(obj_merged, group.by.vars = c("Method",""))

# UMAP after Harmony
obj_merged <- RunUMAP(obj_merged, reduction = "harmony", dims = 1:20)
obj_merged <- FindNeighbors(obj_merged, reduction = "harmony", dims = 1:20)
obj_merged <- FindClusters(obj_merged, resolution = 0.5)

# -----------------------------------------------
# STEP 4: Visualization
# -----------------------------------------------
p_after_method <- DimPlot(obj_merged, reduction = "umap", group.by = "Method") + 
  ggtitle("UMAP After Harmony Integration")

p_after_clusters <- DimPlot(obj_merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("UMAP by Clusters (Post-Integration)")

# Show all UMAPs
p_before / (p_after_method + p_after_clusters)

# Optional: Cluster composition table
table(obj_merged$seurat_clusters, obj_merged$Method)