# Load libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
#install.packages('BiocManager')
#BiocManager::install('glmGamPoi')
# Setup
options(timeout = 1000)
#InstallData("pbmcsca")
options(future.globals.maxSize = 1e9)

# -----------------------------------------------
# STEP 1: Load and Select Methods
# -----------------------------------------------
obj <- LoadData("pbmcsca")
obj <- subset(obj, nFeature_RNA > 1000)

methods_to_use <- c("Smart-seq2", "10x Chromium", "CEL-Seq2")
obj_subset <- subset(obj, subset = Method %in% methods_to_use)
obj_list <- SplitObject(obj_subset, split.by = "Method")

# -----------------------------------------------
# STEP 2: Safe SCTransform and remove NULLs
# -----------------------------------------------
obj_list_sct <- list()

for (method in names(obj_list)) {
  cat("Processing:", method, "\n")
  obj <- obj_list[[method]]
  if (ncol(obj) >= 50) {  # ensure enough cells
    obj <- tryCatch(
      SCTransform(obj, verbose = FALSE),
      error = function(e) { message("Failed SCTransform on ", method); return(NULL) }
    )
    if (!is.null(obj)) {
      obj_list_sct[[method]] <- obj
    }
  } else {
    message("Too few cells in ", method)
  }
}

# -----------------------------------------------
# STEP 3: PCA/UMAP BEFORE INTEGRATION
# -----------------------------------------------
# Merge non-integrated datasets (raw RNA assay)
obj_unintegrated <- merge(obj_list_sct[[1]], y = obj_list_sct[2:length(obj_list_sct)])

# Normalize/scale manually for visualization only
DefaultAssay(obj_unintegrated) <- "RNA"
obj_unintegrated <- NormalizeData(obj_unintegrated)
obj_unintegrated <- FindVariableFeatures(obj_unintegrated)
obj_unintegrated <- ScaleData(obj_unintegrated)
obj_unintegrated <- RunPCA(obj_unintegrated)
obj_unintegrated <- RunUMAP(obj_unintegrated, dims = 1:20)

p_before <- DimPlot(obj_unintegrated, reduction = "umap", group.by = "Method") +
  ggtitle("UMAP Before CCA Integration")

# -----------------------------------------------
# STEP 4: Integration with SCTransform + CCA
# -----------------------------------------------
features <- SelectIntegrationFeatures(object.list = obj_list_sct, nfeatures = 3000)
obj_list_sct <- PrepSCTIntegration(object.list = obj_list_sct, anchor.features = features)

anchors <- FindIntegrationAnchors(
  object.list = obj_list_sct,
  normalization.method = "SCT",
  anchor.features = features
)

obj_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# -----------------------------------------------
# STEP 5: PCA/UMAP/Clusters on Integrated
# -----------------------------------------------
obj_integrated <- RunPCA(obj_integrated)
obj_integrated <- RunUMAP(obj_integrated, dims = 1:20)
obj_integrated <- FindNeighbors(obj_integrated, dims = 1:20)
obj_integrated <- FindClusters(obj_integrated, resolution = 0.5)

p_after_method <- DimPlot(obj_integrated, reduction = "umap", group.by = "Method") +
  ggtitle("UMAP After Integration")
p_after_cluster <- DimPlot(obj_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters After Integration")

# -----------------------------------------------
# STEP 6: Combine Plots
# -----------------------------------------------
p_before / (p_after_method + p_after_cluster)

# Cluster composition table
table(obj_integrated$seurat_clusters, obj_integrated$Method)

