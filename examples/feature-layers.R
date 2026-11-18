# Run after installing Seurcery 0.90 and Seurat:
# Rscript examples/feature-layers.R feature-layers.pdf
# By default, the PDF is saved in the current directory.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Seurcery)
})

data("pbmc_small", package = "SeuratObject")
rds <- pbmc_small
rds$seurat_clusters <- Idents(rds)
gl <- c("LYZ", "MS4A1", "CD3E", "NKG7")
pal <- setNames(c("#3FA55E", "#50B5D8", "#EA8350"), levels(rds$seurat_clusters))

plots <- list(
  categorical = FeatureLayerPlot(
    rds, gl, group.by = "seurat_clusters", reduction = "pca", assay = "RNA"
  ),
  continuous = FeatureLayerPlot(
    rds, gl, group.by = "nCount_RNA", reduction = "pca", assay = "RNA"
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
