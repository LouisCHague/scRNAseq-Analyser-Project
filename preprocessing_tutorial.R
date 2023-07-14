# PREPARE CELL DATA FOR scRNA ANALYSER
# Instructions based on Seurat tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = r"(C:\Users\...\filtered_gene_bc_matrices\hg19)")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Find mitochondrial genes (High % = Might indicate cell death)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Normalizes the data 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# We want to look at how the cell types differ in gene expression
# Most cells express genes at 0 
# Only want to compare genes that are highly variable in their expression
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


# Performs a linear transform (scales) the data
# Shifts the expression of each gene so the mean is 0 and the variance is 1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Performs PCA on the data (Dimensional reduction)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Clusters the cells based on similar gene expression
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# If you know the cell types, you can make cluster ID's 
# If you don't you can use the app's heatmap feature to assist you in finding
# markers that define clusters. 

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$celltype <- Idents(pbmc)

saveRDS(pbmc, file = r"(C:\Users\...)")
