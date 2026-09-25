# Run after installing Seurcery 0.8.7 and Seurat:
# Rscript examples/stacked-gradients.R /path/to/stacked-gradients.pdf
# Without a path, output is written to R's temporary directory.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Seurcery)
})

set.seed(42)
data("pbmc_small", package = "SeuratObject")
rds <- pbmc_small
rds$seurat_clusters <- Idents(rds)
gl <- rownames(rds)[seq_len(15L)]
stopifnot(length(gl) == 15L, !anyDuplicated(gl))
gene_cols <- setNames(
  c(rep("#3FA55E", 2), rep("#50B5D8", 5), rep("#EA8350", 8)),
  gl
)

plots <- list()
for (flip in c(TRUE, FALSE)) {
  p <- Seurat::VlnPlot(
    rds,
    features = gl,
    group.by = "seurat_clusters",
    pt.size = 0,
    stack = TRUE,
    flip = flip,
    add.noise = FALSE
  ) + scale_fill_manual(values = gene_cols)

  prefix <- paste0("flip=", flip)
  plots[[paste(prefix, "violin")]] <-
    p + geom_vln_gradient(alpha_min = 0.05, alpha_max = 1)

  for (guide in c(FALSE, TRUE)) {
    plots[[paste(prefix, "violin then dots; guide=", guide)]] <-
      p +
      geom_vln_gradient(alpha_min = 0.05, alpha_max = 1) +
      geom_dot_gradient(guide = guide)
    plots[[paste(prefix, "dots then violin; guide=", guide)]] <-
      p +
      geom_dot_gradient(guide = guide) +
      geom_vln_gradient(alpha_min = 0.05, alpha_max = 1)
  }
}

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else {
  file.path(tempdir(), "Seurcery-v0.8.7-stacked-gradients.pdf")
}
grDevices::pdf(output, width = 12, height = 11, onefile = TRUE)
tryCatch({
  for (name in names(plots)) {
    message("Rendering: ", name)
    print(plots[[name]])
  }
}, finally = grDevices::dev.off())
message("Saved: ", normalizePath(output))
