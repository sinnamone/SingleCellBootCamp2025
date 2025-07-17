# ---------------------------------------------
# Load libraries
# ---------------------------------------------
library(Seurat)        # Core single-cell analysis
library(SeuratData)    # Seurat example datasets
library(dplyr)         # Data manipulation
library(ggplot2)       # Plotting
library(slingshot)     # Lineage inference
library(SCORPIUS)      # Pseudotime & trajectory analysis
library(viridis)       # Nice color palettes

# ---------------------------------------------
# Load and preprocess PBMC data
# ---------------------------------------------
InstallData("pbmc3k")
pbmc <- LoadData("pbmc3k")

# Quality control
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# SCTransform normalization and basic clustering
pbmc <- SCTransform(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)
Idents(pbmc) <- "seurat_clusters"

# ---------------------------------------------
# Slingshot Trajectory Inference on all cells
# ---------------------------------------------
DimPlot(pbmc,group.by="seurat_clusters")
umap <- Embeddings(pbmc, "umap")      # Extract UMAP coordinates
clust <- Idents(pbmc)                 # Get cluster identities

# Run Slingshot using UMAP and clusters
sling <- slingshot(umap, clusterLabels = clust, reducedDim = "UMAP")
pbmc$pseudotime_slingshot <- slingPseudotime(sling)[,1]  # Save pseudotime

# Visualize Slingshot pseudotime
FeaturePlot(pbmc, features = "pseudotime_slingshot") + ggtitle("Slingshot Pseudotime")


# UMAP plot with Slingshot trajectory overlay
plot(umap, col = viridis(100)[cut(pbmc$pseudotime_slingshot, 100)], pch = 16, asp = 1,
     main = "Slingshot Trajectory (UMAP)")
lines(SlingshotDataSet(sling), lwd = 2, col = 'black')


pbmc$pseudotime_slingshot <- slingPseudotime(sling)[,3]
# UMAP plot with Slingshot trajectory overlay
tr3<- plot(umap, col = viridis(100)[cut(pbmc$pseudotime_slingshot, 100)], pch = 16, asp = 1,
     main = "Slingshot Trajectory (UMAP)") +lines(SlingshotDataSet(sling), lwd = 2, col = 'black')

# ---------------------------------------------
# SCORPIUS Trajectory Inference on all cells
# ---------------------------------------------
# Switch back to RNA assay for SCORPIUS (uses log-normalized data)
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 1000)
hvg <- VariableFeatures(pbmc)

# Extract expression matrix (cells x HVGs)
expr_mat <- t(as.matrix(GetAssayData(pbmc, assay = "RNA", slot = "data")[hvg, ]))

# Dimensionality reduction and trajectory inference
space <- reduce_dimensionality(expr_mat, ndim = 2)
traj <- infer_trajectory(space)

# Plot inferred trajectory
draw_trajectory_plot(space, traj$path, progression_group = traj$time)

# Store pseudotime
pbmc$pseudotime_scorpius <- traj$time

# UMAP visualization
FeaturePlot(pbmc, features = "pseudotime_scorpius") + ggtitle("SCORPIUS Pseudotime")

# ---------------------------------------------
# Subset only clusters 0, 1, and 4 for focused analysis
# ---------------------------------------------
subset_clusters <- c("0", "4", "2")
pbmc_subset <- subset(pbmc, idents = subset_clusters)

# ---------------------------------------------
# Slingshot on subset
# ---------------------------------------------
umap_subset <- Embeddings(pbmc_subset, "umap")
clust_subset <- Idents(pbmc_subset)

sling_subset <- slingshot(umap_subset, clusterLabels = clust_subset, reducedDim = "UMAP")
pbmc_subset$pseudotime_slingshot <- slingPseudotime(sling_subset)[,1]

# Visualize Slingshot pseudotime (subset)
FeaturePlot(pbmc_subset, features = "pseudotime_slingshot") + ggtitle("Subset: Slingshot Pseudotime")
plot(umap_subset, col = viridis(100)[cut(pbmc_subset$pseudotime_slingshot, 100)], pch = 16,
     asp = 1, main = "Slingshot Trajectory (Subset)")
lines(SlingshotDataSet(sling_subset), lwd = 2, col = 'black')

# ---------------------------------------------
# SCORPIUS on subset
# ---------------------------------------------
DefaultAssay(pbmc_subset) <- "RNA"
pbmc_subset <- NormalizeData(pbmc_subset)
pbmc_subset <- FindVariableFeatures(pbmc_subset, nfeatures = 1000)
hvg_subset <- VariableFeatures(pbmc_subset)

# Prepare matrix for SCORPIUS
expr_mat_subset <- t(as.matrix(GetAssayData(pbmc_subset, assay = "RNA", slot = "data")[hvg_subset, ]))
space_subset <- reduce_dimensionality(expr_mat_subset, ndim = 2)
traj_subset <- infer_trajectory(space_subset)

# Visualize SCORPIUS trajectory (subset)
draw_trajectory_plot(space_subset, traj_subset$path, progression_group = traj_subset$time)

# Store pseudotime
pbmc_subset$pseudotime_scorpius <- traj_subset$time
FeaturePlot(pbmc_subset, features = "pseudotime_scorpius") + ggtitle("Subset: SCORPIUS Pseudotime")

# ---------------------------------------------
# Identify trajectory-informative genes (subset)
# ---------------------------------------------
gimp <- gene_importances(
  expr_mat_subset, 
  traj_subset$time, 
  num_permutations = 10, 
  num_threads = 8, 
  ntree = 10000,
  ntree_perm = 1000
) 

# Adjust p-values and select significant genes
gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
gene_sel <- gimp$gene[gimp$qvalue < .05]
expr_sel <- scale_quantile(expr_mat_subset[,gene_sel])  # Quantile normalize

# Plot expression heatmap across pseudotime
draw_trajectory_heatmap(expr_sel, traj_subset$time, pbmc_subset$seurat_clusters)

# ---------------------------------------------
# Extract co-expression gene modules
# ---------------------------------------------
modules <- extract_modules(scale_quantile(expr_sel), traj_subset$time, verbose = F)

# Plot heatmap with modules
draw_trajectory_heatmap(expr_sel, traj_subset$time, pbmc_subset$seurat_clusters, modules,show_labels_row  =T)