# Seurat and DE analysis
library(Seurat)
library(dplyr)
library(SeuratData)
# Enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)     # For human genes
library(enrichplot)
library(msigdbr)          # Optional: MSigDB collections
library(ggplot2)
#install.packages("msigdbr")

# ---------------------------------------------
# STEP 1: Load and preprocess PBMC data
# ---------------------------------------------
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
Idents(pbmc) <- "seurat_clusters"

# ---------------------------------------------
# STEP 2: Marker detection and gene conversion
# ---------------------------------------------
# Find markers for each cluster
all_markers <- FindAllMarkers(pbmc,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

# Background gene list for enrichment
background_genes <- unique(all_markers$gene)
background_entrez <- bitr(background_genes,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

# Load MSigDB Hallmark gene sets
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# Create output directory
outdir <- "/Users/simone/Desktop/pbmc_enrichment"
dir.create(outdir, showWarnings = FALSE)

# ---------------------------------------------
# STEP 3: Enrichment analysis per cluster
# ---------------------------------------------
clusters <- sort(unique(all_markers$cluster))

for (cl in clusters) {
  cat("Processing cluster", cl, "\n")
  
  # Get top marker genes
  top_genes <- all_markers %>%
    filter(cluster == cl) %>%
    top_n(50, wt = avg_log2FC) %>%
    pull(gene)
  
  # Convert to Entrez IDs
  gene_entrez <- bitr(top_genes,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
  
  # Skip if not enough mapped genes
  if (nrow(gene_entrez) < 5) {
    cat("Skipping cluster", cl, "- too few Entrez IDs\n")
    next
  }
  
  # ---------------------
  # GO Biological Process
  ego <- enrichGO(gene = gene_entrez$ENTREZID,
                  universe = background_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  if (!is.null(ego) && nrow(ego) > 0) {
    ego_df <- as.data.frame(ego)
    write.csv(ego_df, file = file.path(outdir, paste0("cluster_", cl, "_GO.csv")), row.names = FALSE)
  }
  
  # ---------------------
  # MSigDB Hallmark enrichment
  emsig <- enricher(gene = gene_entrez$ENTREZID,
                    TERM2GENE = msig_h,
                    universe = background_entrez$ENTREZID)
  
  if (!is.null(emsig) && nrow(emsig) > 0) {
    emsig_df <- as.data.frame(emsig)
    write.csv(emsig_df, file = file.path(outdir, paste0("cluster_", cl, "_MSigDB_Hallmark.csv")), row.names = FALSE)
  }
}
