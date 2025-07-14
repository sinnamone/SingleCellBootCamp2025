# ----------------------------
# STEP 1: Load filtered object
# ----------------------------

# If starting fresh, load the filtered object from Part 2
# pbmc_filtered <- readRDS("pbmc_filtered_seurat.rds")

# ----------------------------
# STEP 2: Normalize the data
# ----------------------------

pbmc_filtered <- NormalizeData(pbmc_filtered)

# ----------------------------
# STEP 3: Identify variable features (HVGs)
# ----------------------------

pbmc_filtered <- FindVariableFeatures(pbmc_filtered, selection.method = "vst", nfeatures = 2000)

# Plot variable features (optional for teaching)
VariableFeaturePlot(pbmc_filtered)

# ----------------------------
# STEP 4: Scale the data
# ----------------------------

# Regress out mitochondrial content and total UMI count (optional but common)
pbmc_filtered <- ScaleData(pbmc_filtered, vars.to.regress = c("nCount_RNA", "percent.mt"))

# ----------------------------
# STEP 5: Run PCA
# ----------------------------

pbmc_filtered <- RunPCA(pbmc_filtered, features = VariableFeatures(pbmc_filtered))

# Examine and visualize PCA results
print(pbmc_filtered[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc_filtered, dims = 1, reduction = "pca")
DimPlot(pbmc_filtered, reduction = "pca")
ElbowPlot(pbmc_filtered)  # use to choose how many PCs to use

# ----------------------------
# STEP 6: Find Neighbors and Clusters
# ----------------------------

# Choose number of PCs (e.g., first 10)
pbmc_filtered <- FindNeighbors(pbmc_filtered, dims = 1:10)

# Find clusters (resolution controls granularity)
pbmc_filtered <- FindClusters(pbmc_filtered, resolution = 0.5)

# View cluster assignments
head(Idents(pbmc_filtered))
table(Idents(pbmc_filtered))

# ----------------------------
# STEP 7: Save object with clustering
# ----------------------------

saveRDS(pbmc_filtered, file = "pbmc_clustered_seurat.rds")


# ----------------------------
# STEP 1: Run UMAP for visualization
# ----------------------------

pbmc_filtered <- RunUMAP(pbmc_filtered, dims = 1:10)

# Plot UMAP colored by cluster
DimPlot(pbmc_filtered, reduction = "umap", label = TRUE, pt.size = 0.5)

# Optionally, plot UMAP colored by metadata (e.g., group)
DimPlot(pbmc_filtered, reduction = "umap", group.by = "group")

# ----------------------------
# STEP 2: Explore individual genes
# ----------------------------

# Feature plot of selected genes (e.g., canonical markers)
FeaturePlot(pbmc_filtered, features = c("CD3D", "MS4A1", "LYZ", "CD8A"))

# ----------------------------
# STEP 3: Find cluster marker genes
# ----------------------------

# Identify markers for all clusters compared to all remaining cells
pbmc_markers <- FindAllMarkers(pbmc_filtered,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

# View top markers per cluster
pbmc_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# ----------------------------
# STEP 4: Save results
# ----------------------------

saveRDS(pbmc_filtered, file = "pbmc_umap_seurat.rds")
write.csv(pbmc_markers, file = "pbmc_cluster_markers.csv", row.names = FALSE)
