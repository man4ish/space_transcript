# Load necessary libraries
library(Seurat)
library(SpatialExperiment)
library(ggplot2)

# Define file paths
data_dir <- "/path/to/spaceranger/output"  # Directory containing the spaceranger output
sample_name <- "sample1"  # Sample ID used in spaceranger count

# Load the spatial transcriptomics data
# In this example, we assume that spaceranger output is in the form of a Seurat-compatible format.
seurat_obj <- Read10X(data.dir = data_dir)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_name)

# Preprocess the data
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

# Perform dimensionality reduction
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Plot UMAP
DimPlot(seurat_obj, reduction = "umap") + ggtitle("UMAP of Spatial Transcriptomics Data")

# Identify clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Plot clusters
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") + ggtitle("Clusters on UMAP")

# Perform spatial analysis if spatial coordinates are available
# Assuming spatial coordinates are included in the data and are accessible
if (!is.null(seurat_obj@images$`Spatial`)) {
  SpatialDimPlot(seurat_obj) + ggtitle("Spatial Distribution")
}

# Save the Seurat object and results
saveRDS(seurat_obj, file = "seurat_analysis_results.rds")

# Optionally, save UMAP and cluster plots
ggsave("umap_plot.png")
ggsave("clusters_plot.png")

