library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
#install.packages("clustree")
library(clustree)
library(dplyr)
library(tibble)
# install.packages('BiocManager')
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)  

# ---------------------------------------------
# STEP 1: Load Data & QC
# ---------------------------------------------
options(timeout = 1000)
#InstallData("pbmc3k")
pbmc <- LoadData("pbmc3k")

# Add mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC Violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ---------------------------------------------
# STEP 2: SCTransform, PCA, UMAP
# ---------------------------------------------
pbmc <- SCTransform(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc)
ElbowPlot(pbmc)  # Choose dimensions (e.g., 1:20)

pbmc <- RunUMAP(pbmc, dims = 1:20)

# ---------------------------------------------
# STEP 3: Clustering at Multiple Resolutions
# ---------------------------------------------
pbmc <- FindNeighbors(pbmc, dims = 1:20)

# Run clustering at multiple resolutions
pbmc <- FindClusters(pbmc, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0))

# Clustree visualization
clustree(pbmc, prefix = "SCT_snn_res.") + ggtitle("Clustree - Cluster Stability")

# Choose a resolution (e.g., 0.4)
Idents(pbmc) <- "SCT_snn_res.0.4"

# UMAP with final clustering
DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.4") +
  ggtitle("UMAP with Clustering (res = 0.4)")

# ---------------------------------------------
# STEP 4: Differential Expression Between Clusters
# ---------------------------------------------
# Compare cluster 0 vs cluster 1 (you can change)
markers_0_vs_1 <- FindMarkers(pbmc, ident.1 = 0, ident.2 = 1, 
                              logfc.threshold = 0.25, min.pct = 0.1)

# Top markers table
head(markers_0_vs_1)

# Add gene names to rownames for clarity
markers_0_vs_1 <- rownames_to_column(markers_0_vs_1, var = "gene")

# ---------------------------------------------
# STEP 5: DEG Visualization
# ---------------------------------------------

# Volcano plot (optional, needs EnhancedVolcano)
EnhancedVolcano(markers_0_vs_1,
                lab = markers_0_vs_1$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Cluster 0 vs Cluster 1',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.5,
                labSize = 4)

# Heatmap of top markers for all clusters
top_markers_all <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- top_markers_all %>% group_by(cluster) %>% top_n(10, avg_log2FC)

DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Feature plot for one or more genes
FeaturePlot(pbmc, features = c("IL7R", "CD3D", "MS4A1"))  # Example markers