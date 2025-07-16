# -------------------------------------------
# Load Libraries
# -------------------------------------------
library(Seurat)
library(SeuratData)
library(dplyr)
library(ggplot2)
library(patchwork)
library(plyr)
# BiocManager::install("SingleR")
library(SingleR)
# BiocManager::install("celldex")
library(celldex)
library(SummarizedExperiment)
# BiocManager::install("scmap")
library(scmap)
# install.packages('SCINA')
library(SCINA)
library(SeuratDisk)

# -------------------------------------------
# Load 10X PBMC data and create Seurat object
# -------------------------------------------
pbmc.data <- Read10X(data.dir = "/Users/simone/Downloads/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "PBMC3K", min.cells = 3, min.features = 200)

# -------------------------------------------
# QC, Normalization, Dimensionality Reduction
# -------------------------------------------
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- SCTransform(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.4)
# Make sure identities are set to clusters
Idents(pbmc) <- "seurat_clusters"

# Subset the object, excluding clusters "7" and "8"
pbmc <- subset(pbmc, idents = setdiff(levels(pbmc), c("7", "8")))

# Reset factor levels (optional but clean)
pbmc$seurat_clusters <- droplevels(pbmc$seurat_clusters)
Idents(pbmc) <- "seurat_clusters"

# -------------------------------------------
# Manual Annotation via Known Marker Genes
# -------------------------------------------
FeaturePlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "NKG7", "PPBP"))

# Define manual annotations
cluster_annotation <- c(
  "0" = "CD4+ T",
  "1" = "CD14+ Monocytes",
  "2" = "B cells",
  "3" = "CD8+ T",
  "4" = "NK cells",
  "5" = "FCGR3A+ Monocytes",
  "6" = "Platelets"
)

pbmc$celltype_manual <- plyr::mapvalues(
  x = as.character(pbmc$seurat_clusters),
  from = names(cluster_annotation),
  to = cluster_annotation
)

# -------------------------------------------
# Save Seurat Object for Azimuth
# -------------------------------------------
# Remove problematic metadata types before saving
pbmc@meta.data <- pbmc@meta.data[, sapply(pbmc@meta.data, function(x) is.atomic(x) || is.factor(x))]

#SaveH5Seurat(pbmc, filename = "~/Desktop/pbmc_for_azimuth.h5seurat", overwrite = TRUE)

# -------------------------------------------
# Automated Annotation with SingleR
# -------------------------------------------
ref <- celldex::HumanPrimaryCellAtlasData()
singleR_result <- SingleR(
  test = GetAssayData(pbmc, slot = "data", assay = "SCT"),
  ref = ref,
  labels = ref$label.main
)
pbmc$celltype_singleR <- singleR_result$labels

# -------------------------------------------
# Automated Annotation with scType (SCINA)
# -------------------------------------------
gene_markers <- list(
  "B cells" = c("MS4A1"),
  "CD4+ T cells" = c("CD3D", "IL7R"),
  "CD8+ T cells" = c("CD8A"),
  "NK cells" = c("NKG7", "GNLY"),
  "Monocytes" = c("CD14", "LYZ"),
  "Dendritic cells" = c("FCER1A"),
  "Platelets" = c("PPBP")
)

expr <- as.matrix(GetAssayData(pbmc, assay = "SCT", slot = "data"))
scina_result <- SCINA(expr, signatures = gene_markers, max_iter = 100)
pbmc$celltype_scType <- scina_result$cell_labels

# -------------------------------------------
# Automated Annotation with scmap
# -------------------------------------------
sce <- as.SingleCellExperiment(pbmc)
colData(sce)$cluster <- as.character(pbmc$seurat_clusters)
rowData(sce)$feature_symbol <- rownames(sce)

# Build index on reference (self-reference for demo)
sce <- selectFeatures(sce, n_features = 500, suppress_plot = TRUE)
sce <- indexCluster(sce, cluster_col = "cluster")

# Project and predict
sce_test <- selectFeatures(sce, n_features = 500, suppress_plot = TRUE)
sce_test <- scmapCluster(
  projection = sce_test,
  index_list = list(ref = metadata(sce)$scmap_cluster_index)
)

pbmc$celltype_scmap <- plyr::mapvalues(
  x = as.character(sce_test$scmap_cluster_labs),
  from = names(cluster_annotation),
  to = cluster_annotation
)

# -------------------------------------------
# Final Overview Plot
# -------------------------------------------
p1 <- DimPlot(pbmc, group.by = "celltype_manual", label = TRUE) + ggtitle("Manual Annotation")
p2 <- DimPlot(pbmc, group.by = "celltype_singleR", label = TRUE) + ggtitle("SingleR")
p3 <- DimPlot(pbmc, group.by = "celltype_scType", label = TRUE) + ggtitle("scType (SCINA)")
p4 <- DimPlot(pbmc, group.by = "celltype_scmap", label = TRUE) + ggtitle("scmap")

(p1 | p2) / (p3 | p4)


counts <- as.matrix(GetAssayData(pbmc, assay = "RNA", slot = "counts"))
write.csv(counts, "~/Desktop/pbmc_counts.csv",quote = F)


#conda create -n celltypist310 python=3.10 -y        
#conda activate celltypist310  
#celltypist -i ~/Desktop/pbmc_counts.csv --model Immune_All_High.pkl --outdir ~/Desktop/   


pbmc$celltype_celltypist <- read.csv("~/Desktop/predicted_labels.csv")$predicted_labels
p5 <- DimPlot(pbmc, group.by = "celltype_celltypist", label = TRUE) + ggtitle("Celltypyst")
