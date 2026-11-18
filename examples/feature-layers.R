# Run after installing Seurcery, Seurat and the PBMC3K data package (see README).
# Rscript examples/feature-layers.R feature-layers.pdf
# By default, the PDF is saved in the current directory.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Seurcery)
})

data("pbmc3k.final", package = "pbmc3k.SeuratData")
rds <- SeuratObject::UpdateSeuratObject(pbmc3k.final)
gl <- c("CD3D", "NKG7", "MS4A1", "LYZ")
clusters <- levels(rds$seurat_clusters)
pal <- setNames(grDevices::hcl.colors(length(clusters), "Dark 3"), clusters)

plots <- list(
  umap = FeatureLayerPlot(
    rds, gl, group.by = "seurat_clusters", reduction = "umap", assay = "RNA"
  ),
  pca = FeatureLayerPlot(
    rds, gl, group.by = "seurat_clusters", reduction = "pca", assay = "RNA"
  ),
  manual_color = FeatureLayerPlot(
    rds, gl, group.by = "seurat_clusters", reduction = "pca",
    group.scale = scale_color_manual(values = pal)
  ),
  manual_fill = FeatureLayerPlot(
    rds, gl, group.by = "seurat_clusters", reduction = "pca",
    group.scale = scale_fill_manual(values = pal)
  ),
  compact = FeatureLayerPlot(
    rds, gl, group.by = "seurat_clusters", reduction = "pca",
    order = TRUE, gap.x = 0.06, gap.y = 0.11, title.gap = 0.04
  )
)

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else "feature-layers.pdf"
grDevices::pdf(output, width = 180 / 25.4, height = 101 / 25.4, onefile = TRUE)
tryCatch({
  for (name in names(plots)) {
    message("Rendering: ", name)
    print(plots[[name]])
  }
}, finally = grDevices::dev.off())
message("Saved: ", normalizePath(output))
