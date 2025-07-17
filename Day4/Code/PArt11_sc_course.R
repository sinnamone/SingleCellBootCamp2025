# Load required libraries
suppressMessages({
  library(Seurat)                     # For single-cell RNA-seq analysis
  library(Signac)                    # For scATAC-seq analysis in Seurat
  library(EnsDb.Hsapiens.v86)        # Gene annotation database
  library(BSgenome.Hsapiens.UCSC.hg38)  # Human genome reference
  library(ggplot2)                   # Visualization
  library(dplyr)                     # Data manipulation
})

# Load 10X Multiome filtered matrix (contains both RNA and ATAC)
inputdata.10x <- Read10X("/Users/simone/Downloads/filtered_feature_bc_matrix/")

# Extract RNA and ATAC separately
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Keep only valid peaks that are real genomic intervals
atac_counts <- atac_counts[grep("^chr", rownames(atac_counts)), ]

# Create Seurat object for RNA modality
obj.rna <- CreateSeuratObject(counts = rna_counts, assay = "RNA")

# Calculate mitochondrial percentage
obj.rna[["percent.mt"]] <- PercentageFeatureSet(obj.rna, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(obj.rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Filter low-quality cells
obj.rna <- subset(obj.rna, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

# Define ATAC fragments path
fragpath <- "/Users/simone/Downloads/pbmc_unsorted_3k_atac_fragments.tsv.gz"

# Create ATAC chromatin assay
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  min.cells = 1,
  genome = "hg38",
  fragments = fragpath
)

# Create ATAC Seurat object
obj.atac <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")

# Load gene annotations and convert to UCSC style
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'

# Add gene annotations to ATAC object
Annotation(obj.atac) <- annotations

# Subset only shared cells between RNA and ATAC
common.cells <- intersect(colnames(obj.rna), colnames(obj.atac))

obj.rna <- subset(obj.rna, cells = common.cells)
obj.atac <- subset(obj.atac, cells = common.cells)

# RNA: Normalization, PCA, UMAP
obj.rna <- obj.rna |>
  NormalizeData(verbose=FALSE) |>
  FindVariableFeatures(nfeatures=3000, verbose=FALSE) |>
  ScaleData(verbose=FALSE) |>
  RunPCA(verbose=FALSE) |>
  RunUMAP(dims=1:30, verbose=FALSE)

# ATAC: TF-IDF, LSI, UMAP
obj.atac <- obj.atac |>
  RunTFIDF() |>
  FindTopFeatures(min.cutoff = 'q0') |>
  RunSVD() |>
  RunUMAP(reduction = "lsi", dims = 2:30)

# RNA UMAP
p1 <- DimPlot(obj.rna, reduction = "umap", group.by = "orig.ident", label = TRUE) + ggtitle("RNA")

# ATAC UMAP
p2 <- DimPlot(obj.atac, reduction = "umap", group.by = "orig.ident", label = TRUE) + ggtitle("ATAC")

# Show side-by-side
p1 + p2

