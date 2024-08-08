# Load necessary libraries
if (!requireNamespace("Seurat", quietly = TRUE)) {
    install.packages("Seurat")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
}

library(Seurat)
library(ggplot2)
library(patchwork)

# Define paths to the Cell Ranger output
# Replace these paths with the actual paths to your Cell Ranger output files
output_dir <- "/path/to/your/cellranger_output"
spatial_data_path <- file.path(output_dir, "filtered_feature_bc_matrix")

# Load the Cell Ranger output
seurat_object <- Read10X(data.dir = spatial_data_path)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = seurat_object, project = "SpatialTranscriptomics")

# Perform initial quality control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj)

# Identify variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# Perform PCA
seurat_obj <- RunPCA(seurat_obj)

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj)

# Run UMAP for dimensionality reduction
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Visualize the clusters
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "ident") + ggtitle("UMAP of Spatial Transcriptomics Data")
p2 <- SpatialFeaturePlot(seurat_obj, features = c("GeneA", "GeneB")) + ggtitle("Spatial Distribution of Gene Expression")

# Combine plots
combined_plot <- p1 + p2

# Save the plots
ggsave("umap_clusters.png", plot = p1)
ggsave("spatial_features.png", plot = p2)
ggsave("combined_plot.png", plot = combined_plot)

# Save the Seurat object
saveRDS(seurat_obj, file = "seurat_object.rds")

# Print completion message
print("Analysis completed. Results saved to files.")

