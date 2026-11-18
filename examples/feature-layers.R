# Run after installing Seurcery 0.90 and Seurat:
# Rscript examples/feature-layers.R /path/to/feature-layers.pdf
# Without a path, output is written to R's temporary directory.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Seurcery)
})

data("pbmc_small", package = "SeuratObject")
rds <- pbmc_small
rds$cell_type <- Idents(rds)
gl <- c("LYZ", "MS4A1", "CD3E", "NKG7")
pal <- setNames(c("#3FA55E", "#50B5D8", "#EA8350"), levels(rds$cell_type))

plots <- list(
  categorical = FeatureLayerPlot(
    rds, gl, group.by = "cell_type", reduction = "pca", assay = "RNA"
  ),
  continuous = FeatureLayerPlot(
    rds, gl, group.by = "nCount_RNA", reduction = "pca", assay = "RNA"
  ),
  manual_color = FeatureLayerPlot(
    rds, gl, group.by = "cell_type", reduction = "pca",
    group.scale = scale_color_manual(values = pal)
  ),
  manual_fill = FeatureLayerPlot(
    rds, gl, group.by = "cell_type", reduction = "pca",
    group.scale = scale_fill_manual(values = pal)
  ),
  compact = FeatureLayerPlot(
    rds, gl, group.by = "cell_type", reduction = "pca",
    order = TRUE, gap.x = 0.06, gap.y = 0.11, title.gap = 0.04
  )
)

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else {
  file.path(tempdir(), "Seurcery-v0.90-feature-layers.pdf")
}
grDevices::pdf(output, width = 180 / 25.4, height = 101 / 25.4, onefile = TRUE)
tryCatch({
  for (name in names(plots)) {
    message("Rendering: ", name)
    print(plots[[name]])
  }
}, finally = grDevices::dev.off())
message("Saved: ", normalizePath(output))
