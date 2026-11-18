layer_test_object <- function() {
  obj <- gradient_pbmc_small()
  obj$layer_group <- factor(rep(c("A", "B"), length.out = ncol(obj)))
  obj
}

layer_points <- function(p) {
  Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)
}

test_that("feature layers use normalized data and the selected embedding", {
  obj <- layer_test_object()
  features <- c("LYZ", "MS4A1")
  values <- as.numeric(t(as.matrix(SeuratObject::LayerData(
    obj[["RNA"]], layer = "data"
  )[features, ])))
  for (reduction in c("pca", "tsne")) {
    xy <- SeuratObject::Embeddings(obj[[reduction]])[, 1:2]
    xy <- sweep(xy, 2, apply(xy, 2, function(v) mean(range(v)))) /
      max(apply(xy, 2, function(v) diff(range(v))))
    p <- FeatureLayerPlot(obj, features, "layer_group", reduction = reduction,
                          order = FALSE, ncol = 1)
    points <- layer_points(p)
    expect_equal(points[[1]]$data$x, unname(xy[, 1]))
    expect_equal(points[[1]]$data$y, unname(xy[, 2]))
    expect_equal(points[[2]]$data$value, values)
    expect_true(any(points[[2]]$data$value == 0))
    expect_equal(nrow(points[[2]]$data), length(features) * ncol(obj))
    expect_s3_class(render_gradient_plot(p), "gtable")

    sorted <- FeatureLayerPlot(obj, features, "layer_group", reduction = reduction)
    expect_equal(layer_points(sorted)[[2]]$data$value, sort(values))
    expect_equal(layer_points(sorted)[[1]]$data, points[[1]]$data)
  }
})

test_that("metadata scales retain Seurat colours and separate expression guides", {
  skip_if_not_installed("Seurat")
  obj <- layer_test_object()
  features <- c("LYZ", "MS4A1")
  p <- FeatureLayerPlot(obj, features, "layer_group", reduction = "pca")
  built <- ggplot2::ggplot_build(p)
  point_ids <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"),
                           logical(1)))
  dim_plot <- Seurat::DimPlot(obj, group.by = "layer_group", reduction = "pca",
                             combine = FALSE)[[1]]
  expect_equal(built$data[[point_ids[1]]]$colour,
               ggplot2::ggplot_build(dim_plot)$data[[1]]$colour)
  pal <- c(A = "#207040", B = "#8040B0")
  for (scale in list(ggplot2::scale_colour_manual(values = pal),
                     ggplot2::scale_fill_manual(values = pal))) {
    manual <- FeatureLayerPlot(obj, features, "layer_group", reduction = "pca",
                               group.scale = scale)
    mb <- ggplot2::ggplot_build(manual)
    aes <- scale$aesthetics
    other <- if (aes == "colour") "fill" else "colour"
    expect_equal(mb$data[[point_ids[1]]][[aes]],
                 unname(pal[as.character(obj$layer_group)]))
    expect_equal(mb$data[[point_ids[2]]][[other]],
                 built$data[[point_ids[2]]]$fill)
  }

  numeric <- FeatureLayerPlot(obj, features, "nCount_RNA", reduction = "pca")
  nb <- ggplot2::ggplot_build(numeric)
  expect_equal(layer_points(numeric)[[1]]$data$group, unname(sort(obj$nCount_RNA)))
  expect_equal(layer_points(numeric)[[1]]$data[, c("x", "y")],
               layer_points(p)[[1]]$data[order(obj$nCount_RNA), c("x", "y")])
  raw <- FeatureLayerPlot(obj, features, "nCount_RNA", reduction = "pca", order = FALSE)
  expect_equal(layer_points(raw)[[1]]$data$group, unname(obj$nCount_RNA))
  expect_equal(layer_points(numeric)[[2]]$data, layer_points(p)[[2]]$data)
  expect_false(nb$plot$scales$get_scales("colour")$is_discrete())
  expect_false(nb$plot$scales$get_scales("fill")$is_discrete())
  colours <- nb$data[[point_ids[1]]]$colour
  expect_equal(toupper(colours[c(1, length(colours))]), c("#FFF5F0", "#67000D"))
  grobs <- walk_grobs(render_gradient_plot(numeric))
  bars <- Filter(function(g) inherits(g, "gtable") && "bar" %in% g$layout$name, grobs)
  expect_length(bars, 2)
  labels <- unlist(lapply(Filter(function(g) inherits(g, "text"), grobs), `[[`, "label"))
  expect_true(all(c("nCount_RNA", "Expression level", "PCA") %in% labels))
})

test_that("spacing and the embedding marker preserve non-overlapping geometry", {
  obj <- layer_test_object()
  features <- c("LYZ", "MS4A1", "CD3E", "NKG7")
  for (ncol in 1:4) {
    p <- FeatureLayerPlot(obj, features, "layer_group", reduction = "pca",
                          ncol = ncol, gap.x = 0.06, gap.y = 0.11)
    frames <- Filter(function(l) inherits(l$geom, "GeomPolygon"), p$layers)[[1]]$data
    boxes <- lapply(split(frames, frames$panel), function(d) {
      c(xmin = min(d$x), xmax = max(d$x), ymin = min(d$y), ymax = max(d$y))
    })
    for (i in seq_len(length(boxes) - 1L)) {
      for (j in seq.int(i + 1L, length(boxes))) {
        a <- boxes[[i]]
        b <- boxes[[j]]
        separate <- a["xmax"] < b["xmin"] || b["xmax"] < a["xmin"] ||
          a["ymax"] + 0.10 <= b["ymin"] || b["ymax"] + 0.10 <= a["ymin"]
        expect_true(unname(separate))
      }
    }
    built <- ggplot2::ggplot_build(p)
    segment_ids <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"),
                               logical(1)))
    marker <- Filter(function(d) nrow(d) == 2L, built$data[segment_ids])[[1]]
    delta <- cbind(marker$xend - marker$x, marker$yend - marker$y)
    expect_equal(sum(delta[1, ] * delta[2, ]), 0)
    expect_equal(sum(delta[1, ]^2), sum(delta[2, ]^2))
    expect_true(all(marker$colour == "black"))
    expect_equal(p$coordinates$ratio, 1)
  }
  expect_error(FeatureLayerPlot(obj, features, "layer_group", gap.x = 0), "gap.x")
  expect_error(FeatureLayerPlot(obj, features, "layer_group", gap.y = 0.09), "gap.y")
  expect_error(FeatureLayerPlot(obj, features, "layer_group", title.gap = 0.03), "title.gap")
  expect_error(FeatureLayerPlot(obj, features, "layer_group", title.gap = 0.2), "gap.y")
})

test_that("absent inputs and invalid coordinates fail before plotting", {
  obj <- layer_test_object()
  expect_error(FeatureLayerPlot(obj, "missing_gene", "layer_group", reduction = "pca"),
               "[Ff]eature|[Gg]ene")
  expect_error(FeatureLayerPlot(obj, "LYZ", "missing_group", reduction = "pca"), "metadata")
  expect_error(FeatureLayerPlot(obj, "LYZ", "layer_group", assay = "missing_assay"), "assay")
  expect_error(FeatureLayerPlot(obj, "LYZ", "layer_group", reduction = "missing_reduction"),
               "reduction")
  bad <- obj
  suppressWarnings(bad[["RNA"]] <- SeuratObject::CreateAssay5Object(
    counts = SeuratObject::LayerData(obj[["RNA"]], layer = "counts")
  ))
  expect_error(FeatureLayerPlot(bad, "LYZ", "layer_group", reduction = "pca"), "data layer")
  for (value in c(NA_real_, NaN, Inf)) {
    bad <- obj
    bad[["pca"]]@cell.embeddings[1, 1] <- value
    expect_error(FeatureLayerPlot(bad, "LYZ", "layer_group", reduction = "pca"),
                 "[Ee]mbedding|[Cc]oordinate")
  }
})

test_that("missing metadata is reported and retained without numbered NA labels", {
  obj <- layer_test_object()
  obj$layer_group[1] <- NA
  expect_warning(p <- FeatureLayerPlot(obj, "LYZ", "layer_group", reduction = "pca"),
                 "1.*[Mm]issing|[Mm]issing.*1|1.*non-finite")
  expect_equal(nrow(layer_points(p)[[1]]$data), ncol(obj))
  expect_equal(sum(is.na(layer_points(p)[[1]]$data$group)), 1)
  built <- ggplot2::ggplot_build(p)
  scale <- built$plot$scales$get_scales("colour")
  expect_equal(tail(as.character(scale$get_labels()), 1), "NA")
  for (value in c(NA_real_, NaN, Inf)) {
    obj$numeric_group <- seq_len(ncol(obj))
    obj$numeric_group[1] <- value
    expect_warning(p <- FeatureLayerPlot(obj, "LYZ", "numeric_group", reduction = "pca"),
                   "1.*[Mm]issing|[Mm]issing.*1|1.*non-finite")
    expect_equal(sum(is.na(layer_points(p)[[1]]$data$group)), 1)
  }
  obj$numeric_group <- NA_real_
  expect_error(FeatureLayerPlot(obj, "LYZ", "numeric_group", reduction = "pca"),
               "[Aa]ll|[Nn]o finite|[Nn]o non-missing")
  obj$layer_group <- NA_character_
  expect_error(FeatureLayerPlot(obj, "LYZ", "layer_group", reduction = "pca"),
               "[Aa]ll|[Nn]o finite|[Nn]o non-missing")
})

test_that("missing expression warns without replacing values by zero", {
  obj <- layer_test_object()
  values <- SeuratObject::LayerData(obj[["RNA"]], layer = "data")
  values["LYZ", 1:3] <- c(NA_real_, NaN, Inf)
  SeuratObject::LayerData(obj[["RNA"]], layer = "data") <- values
  expect_warning(p <- FeatureLayerPlot(obj, c("LYZ", "MS4A1"), "layer_group",
                                       reduction = "pca"),
                 "3.*[Mm]issing|[Mm]issing.*3|3.*non-finite")
  dots <- layer_points(p)[[2]]$data
  expect_equal(nrow(dots), 2L * ncol(obj))
  expect_equal(sum(is.na(dots$value)), 3)
  expect_true(any(dots$value == 0, na.rm = TRUE))
  expect_s3_class(render_gradient_plot(p), "gtable")
  values["LYZ", ] <- NA_real_
  SeuratObject::LayerData(obj[["RNA"]], layer = "data") <- values
  expect_error(FeatureLayerPlot(obj, "LYZ", "layer_group", reduction = "pca"),
               "LYZ")
})

test_that("split data layers preserve cell identities and missing layer features", {
  obj <- layer_test_object()
  features <- c("LYZ", "MS4A1")
  baseline <- FeatureLayerPlot(obj, features, "layer_group", reduction = "pca", order = FALSE)
  suppressWarnings(suppressMessages(obj[["RNA"]] <- split(
    methods::as(obj[["RNA"]], "Assay5"), rep(c("a", "b"), each = ncol(obj) / 2)
  )))
  p <- FeatureLayerPlot(obj, features, "layer_group", reduction = "pca", order = FALSE)
  expect_equal(lapply(layer_points(p), `[[`, "data"),
               lapply(layer_points(baseline), `[[`, "data"))
  values <- SeuratObject::LayerData(obj[["RNA"]], layer = "data.b")
  suppressWarnings(SeuratObject::LayerData(obj[["RNA"]], layer = "data.b") <-
    values[rownames(values) != "LYZ", , drop = FALSE])
  expect_warning(p <- FeatureLayerPlot(obj, features, "layer_group", reduction = "pca",
                                       order = FALSE), "[Mm]issing|non-finite")
  dots <- layer_points(p)[[2]]$data
  expect_equal(nrow(dots), 2L * ncol(obj))
  expect_equal(sum(is.na(dots$value)), ncol(obj) / 2)
  expect_equal(dots$x, layer_points(baseline)[[2]]$data$x)
  expect_equal(dots$y, layer_points(baseline)[[2]]$data$y)
})

test_that("all-zero expression uses the low colour and preserves missing values", {
  obj <- layer_test_object()
  values <- SeuratObject::LayerData(obj[["RNA"]], layer = "data")
  values["LYZ", ] <- 0
  values["LYZ", 1] <- NA_real_
  SeuratObject::LayerData(obj[["RNA"]], layer = "data") <- values
  expect_warning(p <- FeatureLayerPlot(obj, "LYZ", "layer_group", reduction = "pca",
                                       feature.cols = c("#123456", "red")), "1.*missing")
  scale <- ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")
  expect_equal(scale$map(c(0, NA_real_)), c("#123456", "grey50"))
  expect_equal(as.numeric(scale$get_breaks()), 0)
})

test_that("ambiguous normalized layers identify the overlapping gene and cell", {
  obj <- layer_test_object()
  assay <- methods::as(obj[["RNA"]], "Assay5")
  values <- SeuratObject::LayerData(assay, layer = "data")
  duplicate <- values[c("LYZ", "MS4A1"), 1:10, drop = FALSE]
  duplicate[,] <- 5
  SeuratObject::LayerData(assay, layer = "data.alt") <- duplicate
  suppressWarnings(obj[["RNA"]] <- assay)
  expect_error(FeatureLayerPlot(obj, "LYZ", "layer_group", reduction = "pca"),
               paste0("Overlapping.*LYZ.*", colnames(duplicate)[1]))
  expect_no_error(FeatureLayerPlot(obj, "NKG7", "layer_group", reduction = "pca"))
})
