# Seurat and DE analysis
library(Seurat)
library(dplyr)
library(SeuratData)
# Enrichment analysis
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)     # For human genes
library(enrichplot)
#install.packages("msigdbr")
library(msigdbr)          # Optional: MSigDB collections
library(ggplot2)



# ---------------------------------------------
# STEP 1: Load and preprocess Seurat object
# ---------------------------------------------
# Load pbmc3k dataset
InstallData("pbmc3k")
pbmc <- LoadData("pbmc3k")

# Quality control
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize and cluster
pbmc <- SCTransform(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# ---------------------------------------------
# STEP 2: Find cluster markers
# ---------------------------------------------
Idents(pbmc) <- "seurat_clusters"

all_markers <- FindAllMarkers(pbmc,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

# Background gene universe
background_genes <- unique(all_markers$gene)
background_entrez <- bitr(background_genes,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

# ---------------------------------------------
# STEP 3: Enrichment analysis (GO, KEGG)
# ---------------------------------------------
# Choose a cluster to analyze
cluster_of_interest <- 0

top_genes <- all_markers %>%
  filter(cluster == cluster_of_interest) %>%
  top_n(50, wt = avg_log2FC)

genes_entrez <- bitr(top_genes$gene,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

# GO Biological Process
ego <- enrichGO(gene = genes_entrez$ENTREZID,
                universe = background_entrez$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# KEGG
ekegg <- enrichKEGG(gene = genes_entrez$ENTREZID,
                    universe = background_entrez$ENTREZID,
                    organism = 'hsa',
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05)

# ---------------------------------------------
# STEP 4: Visualization
# ---------------------------------------------
# UMAP with cluster labels
DimPlot(pbmc, reduction = "umap", label = TRUE) + ggtitle("UMAP of PBMC Clusters")

# Enrichment plots
dotplot(ego, showCategory = 10, title = "GO BP Enrichment - Cluster 0")
dotplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment - Cluster 0")

# Optional: Barplot
barplot(ego, showCategory = 10)

# ---------------------------------------------
# STEP 5: Optional MSigDB enrichment
# ---------------------------------------------
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")

msig_enrichment <- enricher(gene = genes_entrez$ENTREZID,
                            TERM2GENE = msigdb_hallmark %>% dplyr::select(gs_name, entrez_gene),
                            universe = background_entrez$ENTREZID)

dotplot(msig_enrichment, showCategory = 10, title = "MSigDB Hallmark - Cluster 0")



# Export GO results
if (!is.null(ego) && nrow(ego) > 0) {
  ego_df <- as.data.frame(ego)
  write.csv(ego_df, file = file.path(outdir, paste0("cluster_", cl, "_GO.csv")), row.names = FALSE)
}

# Optional: MSigDB enrichment
emsig <- enricher(gene = gene_entrez$ENTREZID,
                  TERM2GENE = msig_h,
                  universe = background_entrez$ENTREZID)

if (!is.null(emsig) && nrow(emsig) > 0) {
  emsig_df <- as.data.frame(emsig)
  write.csv(emsig_df, file = file.path(outdir, paste0("cluster_", cl, "_MSigDB_Hallmark.csv")), row.names = FALSE)
}

# Optional: KEGG
# ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = 'hsa', universe = background_entrez$ENTREZID)
# if (!is.null(ekegg) && nrow(ekegg) > 0) {
#   write.csv(as.data.frame(ekegg), file = file.path(outdir, paste0("cluster_", cl, "_KEGG.csv")), row.names = FALSE)
# }
}
