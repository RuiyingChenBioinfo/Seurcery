test_that("Seurat stacked violins support both flip settings and manual gene colours", {
  skip_if_not_installed("Seurat")
  obj <- gradient_pbmc_small()
  obj$seurat_clusters <- as.character(Seurat::Idents(obj))
  features <- c("LYZ", "MS4A1", "CD3E")
  cols <- stats::setNames(c("#3FA55E", "#50B5D8", "#EA8350"), features)
  for (flip in c(TRUE, FALSE)) {
    p <- suppressWarnings(Seurat::VlnPlot(
      obj, features = features, group.by = "seurat_clusters",
      pt.size = 0, stack = TRUE, flip = flip
    ))
    out <- suppressMessages(p + geom_vln_gradient(
      bin = 20, alpha_min = 0.1, alpha_max = 1, outline = FALSE,
      min_width_frac = 0
    ) + ggplot2::scale_fill_manual(values = unname(cols)))
    records <- panel_records(out)
    expect_length(records, length(features))
    seen_cols <- character()
    for (panel in records) {
      expect_gt(length(panel), 0)
      fills <- unique(vapply(panel, `[[`, character(1), "rgb"))
      alphas <- vapply(panel, `[[`, numeric(1), "alpha")
      expect_length(fills, 1)
      expect_lt(min(alphas), 0.25)
      expect_gt(max(alphas), 0.8)
      seen_cols <- c(seen_cols, fills)
    }
    expect_setequal(seen_cols, unname(cols))
    combined <- out + geom_dot_gradient(guide = FALSE)
    expect_equal(nrow(dot_summary_data(combined)),
                 length(features) * length(unique(obj$seurat_clusters)))
    expect_s3_class(render_gradient_plot(combined), "gtable")
  }
})

test_that("Seurat split violin geometry retains each half", {
  skip_if_not_installed("Seurat")
  obj <- gradient_pbmc_small()
  obj$seurat_clusters <- as.character(Seurat::Idents(obj))
  obj$gradient_test_split <- rep(c("A", "B"), length.out = ncol(obj))
  p <- suppressMessages(suppressWarnings(Seurat::VlnPlot(
    obj, features = "LYZ", group.by = "seurat_clusters",
    split.by = "gradient_test_split", split.plot = TRUE,
    pt.size = 0, combine = FALSE
  )))[[1]]
  cols <- c(A = "#3FA55E", B = "#EA8350")
  out <- suppressMessages(p + ggplot2::scale_fill_manual(values = cols) +
    geom_vln_gradient(bin = 12, outline = TRUE, min_width_frac = 0))
  expect_s3_class(render_gradient_plot(out), "gtable")
  polygons <- panel_records(out)[[1]]
  expect_setequal(unique(vapply(polygons, `[[`, character(1), "rgb")), unname(cols))
  built <- ggplot2::ggplot_build(out)
  vd <- built$data[[1]]
  centres <- built$layout$coord$transform(vd, built$layout$panel_params[[1]])$x
  centres <- sort(unique(centres))
  for (poly in polygons) {
    # A split polygon has an entire vertical edge on a category centre.
    left_edge <- min(poly$x)
    right_edge <- max(poly$x)
    expect_true(any(abs(centres - left_edge) < 1e-6) ||
                  any(abs(centres - right_edge) < 1e-6))
  }
})
