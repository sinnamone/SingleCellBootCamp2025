# ---------------------------------------------
# Load required libraries
# ---------------------------------------------
suppressMessages({
  library(Seurat)                         # scRNA-seq framework
  library(Signac)                         # scATAC-seq with Seurat
  library(EnsDb.Hsapiens.v86)             # Gene annotation
  library(BSgenome.Hsapiens.UCSC.hg38)    # Genome reference
  library(ggplot2)                        # Plots
  library(dplyr)                          # Data wrangling
})
library(SeuratData)
InstallData("pbmcMultiome")
pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")

# ---------------------------------------------
# Load data: 10x Multiome (RNA + ATAC)
# ---------------------------------------------
# Adjust paths accordingly
data.dir <- "../Downloads/filtered_feature_bc_matrix/"
fragpath <- "../Downloads/pbmc_unsorted_3k_atac_fragments.tsv.gz"

inputdata.10x <- Read10X(data.dir = data.dir)
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Keep only peaks with standard chromosome notation
atac_counts <- atac_counts[grep("^chr", rownames(atac_counts)), ]

# ---------------------------------------------
# Create RNA assay object
# ---------------------------------------------
obj.rna <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
obj.rna[["percent.mt"]] <- PercentageFeatureSet(obj.rna, pattern = "^MT-")

# QC plot
VlnPlot(obj.rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0)

# Filter low-quality cells
obj.rna <- subset(obj.rna, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

# ---------------------------------------------
# Create ATAC assay (ChromatinAssay object)
# ---------------------------------------------
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = fragpath,
  min.cells = 1
)

obj.atac <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")

# ---------------------------------------------
# Add gene annotations to ATAC
# ---------------------------------------------

pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
#annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")
#Annotation(pbmc.atac) <- annotations

# ---------------------------------------------
# Subset to shared cells
seqinfo(pbmc.atac) <- seqinfo(pbmc.rna) 
common.cells <- intersect(colnames(pbmc.rna), colnames(pbmc.atac))
obj.rna <- subset(pbmc.rna, cells = common.cells)
obj.atac <- subset(pbmc.atac, cells = common.cells)

# ---------------------------------------------
# RNA analysis: Normalize, PCA, UMAP
# ---------------------------------------------
obj.rna <- obj.rna |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(nfeatures = 3000, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |>
  RunPCA(verbose = FALSE) |>
  RunUMAP(dims = 1:30, verbose = FALSE)

# ---------------------------------------------
# ATAC analysis: TF-IDF, LSI, UMAP
# ---------------------------------------------
obj.atac <- obj.atac |>
  RunTFIDF() |>
  FindTopFeatures(min.cutoff = 'q0') |>
  RunSVD() |>
  RunUMAP(reduction = "lsi", dims = 2:30)

# ---------------------------------------------
# Plot RNA and ATAC UMAPs side-by-side
# ---------------------------------------------
p1 <- DimPlot(obj.rna, reduction = "umap", group.by = "orig.ident", label = TRUE) + ggtitle("RNA")
p2 <- DimPlot(obj.atac, reduction = "umap", group.by = "orig.ident", label = TRUE) + ggtitle("ATAC")
p1 + p2
