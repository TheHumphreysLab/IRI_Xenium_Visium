library(Seurat)
library(Matrix)
library(ggplot2)

## snRNA-seq data (GSE139107)
load("IRIall.RData")
sc_all <- UpdateSeuratObject(IRI.all)

### Read Scanpy files and convert to Seurat Object (check 0_preprocessing_harmony_integration.py)
counts <- readMM("py2r/IRI_all_count.mtx")
counts <- as(counts, "dgCMatrix")
counts <- t(counts)

genes <- read.csv("py2r/IRI_all_genes.csv")
cells <- read.csv("py2r/IRI_all_cells.csv")
rownames(counts) <- genes$genes
colnames(counts) <- cells$cells

meta <- read.csv("py2r/IRI_all_meta.csv")
rownames(meta) <- meta$X
meta <- meta[,-1]
umap <- read.csv("py2r/IRI_all_umap.csv")
rownames(umap)<-colnames(counts)
umap<- umap[,-1]
names(umap)<-c("UMAP_1", "UMAP_2")
umap <- as.matrix(umap)

spatial <- read.csv("py2r/IRI_all_spatial.csv")
rownames(spatial)<-colnames(counts)
spatial<- spatial[,-1]
names(spatial)<-c("x_coords", "y_coords")
spatial <- as.matrix(spatial)

harmony <- read.csv("py2r/IRI_all_pca_harmony.csv")
rownames(harmony)<-colnames(counts)
harmony<- harmony[,-1]
names(harmony)<-paste("PC_harmony", 1:ncol(harmony), sep="_")
harmony <- as.matrix(harmony)

sp_all <- CreateSeuratObject(counts = counts, meta.data = meta,
											min.cells = 0, min.features = 0)
sp_all <- SCTransform(sp_all)
sp_all <- RunPCA(sp_all)

sp_all@reductions$pca@cell.embeddings <- harmony
sp_all@reductions$harmony <- sp_all@reductions$pca
sp_all@reductions$harmony@cell.embeddings <- harmony
sp_all@reductions$umap <- sp_all@reductions$pca
sp_all@reductions$umap@cell.embeddings <- umap
sp_all@reductions$umap@key <- "UMAP_"

sp_all@reductions$spatial <- sp_all@reductions$pca
sp_all@reductions$spatial@cell.embeddings <- spatial
sp_all@reductions$spatial@key <- "spatial_"

anchor <- FindTransferAnchors(
  reference = sc_all,
  query = sp_all,
  reference.reduction = "pca",
  normalization.method = "SCT",
  dims = 1:30
)

predictions <- TransferData(anchorset = anchor, refdata = sc_all$celltype,
    dims = 1:30)
sp_all <- AddMetaData(sp_all, metadata = predictions)

sp_all <- MapQuery(
  anchorset = anchor,
  query = sp_all,
  reference = sc_all,
  refdata = list(
    celltype = "celltype"),
  reference.reduction = "pca", 
  reduction.model = "umap"
)
