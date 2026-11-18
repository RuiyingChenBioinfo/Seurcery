# Run after installing Seurcery, Seurat and the PBMC3K data package (see README).
# Rscript examples/stacked-gradients.R stacked-gradients.pdf
# By default, the PDF is saved in the current directory.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Seurcery)
})

data("pbmc3k.final", package = "pbmc3k.SeuratData")
rds <- SeuratObject::UpdateSeuratObject(pbmc3k.final)
gl <- c("CD3D", "NKG7", "MS4A1", "LYZ")
gene_cols <- c(CD3D = "#3FA55E", NKG7 = "#50B5D8", MS4A1 = "#EA8350", LYZ = "#8B72BE")

p <- Seurat::VlnPlot(
  rds, features = gl, group.by = "seurat_clusters", pt.size = 0,
  stack = TRUE, flip = TRUE, add.noise = FALSE
) + scale_fill_manual(values = gene_cols)

plots <- list(
  violin = p + geom_vln_gradient(),
  dots = p + geom_dot_gradient(guide = FALSE, y_axis_min = -2),
  combined = p + geom_vln_gradient() + geom_dot_gradient(guide = TRUE, y_axis_min = -2)
)

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else "stacked-gradients.pdf"
grDevices::pdf(output, width = 8, height = 5.5, onefile = TRUE)
tryCatch({
  for (name in names(plots)) {
    message("Rendering: ", name)
    print(plots[[name]])
  }
}, finally = grDevices::dev.off())
message("Saved: ", normalizePath(output))
