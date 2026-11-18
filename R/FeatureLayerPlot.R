#' Plot expression layers on a shared Seurat embedding
#'
#' Places a metadata reference map beside projected gene panels using the same
#' cells and embedding. Expression comes from the selected assay's normalized
#' `data` layer and uses one shared colour scale labelled "Expression level".
#' Continuous metadata uses an embedded Reds palette; categorical metadata uses
#' ggplot2's discrete colours, as in Seurat's default `DimPlot()`.
#'
#' The selected reduction must already exist in the object. Its first two
#' dimensions are used, with an equal-arm embedding marker at the lower left.
#' Missing or non-finite coordinates stop plotting. Partial missing metadata
#' or expression values produce a counted warning and use the scale's
#' missing-value colour (`grey50` by default); they are never replaced by zero.
#' Entirely missing metadata or an entirely missing requested feature is an
#' error. Split Seurat v5 `data.*` layers are read without removing cells whose
#' requested expression values are all missing. If two normalized layers cover
#' the same requested gene and cell, the function reports the ambiguity.
#' When all observed expression is zero, points use the first expression colour.
#'
#' @param object A Seurat object with normalized expression and a stored reduction.
#' @param features Character vector of feature names in the normalized data.
#' @param group.by Name of a metadata column for the reference map. The
#'   `seurat_clusters` column is categorical even when stored as numbers.
#'   Its original IDs appear on the map and in the legend. Other numeric
#'   columns are continuous; factor levels that are integer IDs keep their labels.
#' @param reduction Name of the stored reduction, such as `"umap"` or `"pca"`.
#'   Default is `"umap"`.
#' @param assay Assay containing normalized expression. Defaults to the object's
#'   active assay.
#' @param ncol Number of columns of gene panels. Must be a positive integer.
#' @param pt.size Point size passed to `ggplot2::geom_point()`.
#' @param group.scale Optional ggplot2 colour or fill scale for the reference map.
#'   For categories, use `scale_color_manual()` or `scale_fill_manual()` with
#'   values named by metadata level. Continuous metadata requires a continuous scale.
#' @param feature.cols Colours defining the shared expression gradient.
#' @param text.color Colour of titles, labels, and legend text. Default is black.
#' @param order Logical. Draw larger expression values last. Also applies to
#'   continuous metadata. Default `TRUE` preserves visible high-expression points;
#'   `FALSE` keeps the embedding's cell order.
#' @param gap.x Horizontal gap between gene panels in normalized plot coordinates.
#'   Must be greater than zero. Default is `0.18`.
#' @param gap.y Vertical gap between gene panels in normalized plot coordinates.
#'   Must be at least `max(0.10, title.gap + 0.037)` to reserve title space.
#'   Default is `0.16`.
#' @param title.gap Distance from each gene panel's top to its title centre in
#'   normalized plot coordinates. Must be at least `0.037`. Default is `0.04`.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Install the pbmc3k dataset as shown in the tutorial.
#' data("pbmc3k.final", package = "pbmc3k.SeuratData")
#' rds <- SeuratObject::UpdateSeuratObject(pbmc3k.final)
#' FeatureLayerPlot(
#'   rds, c("CD3D", "NKG7", "MS4A1", "LYZ"),
#'   group.by = "seurat_clusters", reduction = "umap"
#' )
#' }
#' @importFrom rlang .data
#' @export
FeatureLayerPlot <- function(object, features, group.by, reduction = "umap",
                             assay = SeuratObject::DefaultAssay(object),
                             ncol = 2L, pt.size = 0.83, group.scale = NULL,
                             feature.cols = c("lightgrey", "blue"),
                             text.color = "black", order = TRUE,
                             gap.x = 0.18, gap.y = 0.16, title.gap = 0.04) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }
  if (!is.character(features) || !length(features) ||
      anyNA(features) || any(!nzchar(features))) {
    stop("`features` must contain non-missing feature names.", call. = FALSE)
  }
  features <- unique(features)
  for (name in c("group.by", "assay", "reduction")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L ||
        is.na(value) || !nzchar(value)) {
      stop("`", name, "` must be one non-missing name.", call. = FALSE)
    }
  }
  if (!.seurcery_scalar(ncol) || ncol < 1 || ncol != floor(ncol)) {
    stop("`ncol` must be a positive integer.", call. = FALSE)
  }
  if (!.seurcery_scalar(pt.size) || pt.size <= 0) {
    stop("`pt.size` must be a positive finite number.", call. = FALSE)
  }
  if (!is.logical(order) || length(order) != 1L || is.na(order)) {
    stop("`order` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!.seurcery_scalar(gap.x) || gap.x <= 0) {
    stop("`gap.x` must be greater than 0 to keep gene panels separate.", call. = FALSE)
  }
  if (!.seurcery_scalar(title.gap) || title.gap < 0.037) {
    stop("`title.gap` must be at least 0.037 to keep gene titles above panels.", call. = FALSE)
  }
  title.space <- max(0.10, title.gap + 0.037)
  if (!.seurcery_scalar(gap.y) || gap.y < title.space) {
    stop("`gap.y` must be at least ", signif(title.space, 4),
         " to keep gene panels and title space separate.", call. = FALSE)
  }
  meta <- object[[]]
  if (!group.by %in% names(meta)) {
    stop("Unknown metadata column: ", group.by, call. = FALSE)
  }
  if (!assay %in% SeuratObject::Assays(object)) {
    stop("Unknown assay: ", assay, call. = FALSE)
  }
  if (!reduction %in% SeuratObject::Reductions(object)) {
    stop("Unknown reduction: ", reduction, call. = FALSE)
  }
  rna <- object[[assay]]
  layers <- SeuratObject::Layers(rna)
  layers <- layers[grepl("^data($|[.])", layers)]
  if (!length(layers)) {
    stop("The selected assay needs a normalized data layer.", call. = FALSE)
  }
  layer.features <- lapply(layers, function(layer) {
    intersect(features, SeuratObject::Features(rna, layer = layer))
  })
  missing <- setdiff(features, unlist(layer.features, use.names = FALSE))
  if (length(missing)) {
    stop("Missing features in normalized data: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  xy <- SeuratObject::Embeddings(object[[reduction]])
  if (ncol(xy) < 2L) {
    stop("The selected reduction needs at least two dimensions.", call. = FALSE)
  }
  xy <- xy[, 1:2, drop = FALSE]
  cells <- intersect(rownames(xy), colnames(rna))
  if (!length(cells)) {
    stop("No cells are shared by the selected assay and reduction.", call. = FALSE)
  }
  if (length(cells) < nrow(xy)) {
    warning("Omitting ", nrow(xy) - length(cells),
            " embedding cells absent from the selected assay.", call. = FALSE)
  }
  xy <- xy[cells, , drop = FALSE]
  if (any(!is.finite(xy))) {
    stop("Reduction '", reduction, "' contains ", sum(!is.finite(xy)),
         " missing or non-finite coordinates.", call. = FALSE)
  }
  if (max(apply(xy, 2, function(v) diff(range(v)))) == 0) {
    stop("The selected embedding has no coordinate range.", call. = FALSE)
  }
  # Keep missing expression in split layers instead of removing cells.
  if (inherits(rna, "StdAssay")) {
    if (length(layers) > 1L) {
      layer.cells <- lapply(layers, function(layer) {
        intersect(cells, SeuratObject::Cells(rna, layer = layer))
      })
      for (gene in features) {
        has.gene <- vapply(layer.features, function(x) gene %in% x, logical(1))
        ids <- unlist(layer.cells[has.gene], use.names = FALSE)
        dup <- anyDuplicated(ids)
        if (dup) {
          stop("Overlapping normalized layers for feature '", gene,
               "' in cell '", ids[dup], "'.", call. = FALSE)
        }
      }
    }
    expr <- SeuratObject::FetchData(rna, vars = features, cells = cells,
                                    layer = "^data($|[.])", clean = FALSE)
  } else {
    expr <- SeuratObject::FetchData(rna, vars = features, cells = cells, layer = "data")
  }
  expr <- expr[cells, features, drop = FALSE]
  bad.expr <- !is.finite(as.matrix(expr))
  missing <- features[colSums(!bad.expr) == 0L]
  if (length(missing)) {
    stop("All expression values are missing or non-finite for: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (any(bad.expr)) {
    warning("Found ", sum(bad.expr), " missing or non-finite expression values; ",
            "using the scale's missing-value colour.", call. = FALSE)
    expr[bad.expr] <- NA_real_
  }
  if (all(expr == 0, na.rm = TRUE)) feature.cols <- rep(feature.cols[1L], 2L)
  groups <- meta[cells, group.by]
  is.cluster <- group.by == "seurat_clusters"
  continuous <- is.numeric(groups) && !is.cluster
  bad.group <- if (is.numeric(groups)) !is.finite(groups) else is.na(groups)
  if (all(bad.group)) {
    stop("All values are missing or non-finite in metadata '", group.by, "'.",
         call. = FALSE)
  }
  if (any(bad.group)) {
    warning("Found ", sum(bad.group), " missing or non-finite values in metadata '",
            group.by, "'; using the scale's missing-value colour.", call. = FALSE)
    groups[bad.group] <- NA
  }
  if (!continuous) groups <- droplevels(as.factor(groups))
  keys <- levels(groups)
  original.ids <- !continuous && (is.cluster || all(grepl("^-?[0-9]+$", keys)))

  # Keep the original aspect ratio in the reference map.
  xy <- sweep(xy, 2, apply(xy, 2, function(v) mean(range(v)))) /
    max(apply(xy, 2, function(v) diff(range(v))))
  ref <- data.frame(x = xy[, 1], y = xy[, 2], group = groups)
  centres <- data.frame(x = numeric(), y = numeric(), id = integer())
  if (!continuous) {
    centres <- stats::aggregate(cbind(x, y) ~ group, ref, stats::median)
    centres$id <- if (original.ids) as.character(centres$group) else match(centres$group, keys)
    # Put the number beside clusters narrower than the label.
    cluster.width <- tapply(ref$x, ref$group, stats::IQR)
    narrow <- cluster.width[as.character(centres$group)] < 0.025
    centres$x[narrow] <- centres$x[narrow] + 0.043
  }

  # Every expression panel uses the same affine projection.
  project <- function(x, y) data.frame(x = x + 0.22 * y, y = 0.075 * x + 0.49 * y)
  xr <- range(ref$x) + c(-0.04, 0.04)
  yr <- range(ref$y) + c(-0.04, 0.04)
  box <- project(xr[c(1, 2, 2, 1)], yr[c(1, 1, 2, 2)])
  uv <- project(ref$x, ref$y)
  nw <- min(as.integer(ncol), length(features))
  nh <- ceiling(length(features) / nw)
  plate.width <- diff(range(box$x))
  plate.height <- diff(range(box$y))
  pitch.x <- plate.width + gap.x
  pitch.y <- plate.height + gap.y
  gene.left <- max(ref$x) + 0.28
  gene.right <- gene.left + (nw - 1) * pitch.x + plate.width
  first.y <- (nh - 1) * pitch.y / 2
  i <- seq_along(features) - 1L
  dx <- gene.left - min(box$x) + (i %% nw) * pitch.x
  dy <- first.y - (i %/% nw) * pitch.y
  frames <- do.call(rbind, lapply(seq_along(features), function(j) {
    data.frame(x = box$x + dx[j], y = box$y + dy[j], panel = j)
  }))
  dots <- do.call(rbind, lapply(seq_along(features), function(j) {
    data.frame(x = uv$x + dx[j], y = uv$y + dy[j], value = expr[[j]])
  }))
  if (order) {
    dots <- dots[base::order(dots$value, na.last = FALSE), ]
    if (continuous) ref <- ref[base::order(ref$group, na.last = FALSE), ]
  }
  titles <- data.frame(x = min(box$x) + dx, y = max(box$y) + dy + title.gap,
                       label = features)
  title.y <- max(titles$y, max(ref$y) + 0.10)
  titles <- rbind(data.frame(x = min(ref$x) - 0.035, y = title.y,
                             label = group.by), titles)
  bottom <- min(frames$y, ref$y) - 0.025
  left <- min(ref$x) - 0.06
  right <- gene.right + 0.035
  embedding.label <- switch(tolower(reduction), umap = "UMAP", pca = "PCA",
                             tsne = "t-SNE", reduction)
  axis.x <- left + 0.025
  axis.y <- bottom - 0.17
  axis.length <- 0.15

  # Manual scales use the original level names for matching.
  if (is.null(group.scale)) {
    if (continuous) {
      # Matplotlib Reds: RGB interpolation and a 256-colour lookup table.
      reds <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A",
                "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
      rgb <- grDevices::colorRamp(reds, space = "rgb")(seq(0, 1, length.out = 256))
      reds <- grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255)
      group.scale <- ggplot2::scale_colour_gradientn(colours = reds)
      group.scale$palette <- function(x) reds[pmin(256L, floor(x * 256) + 1L)]
    } else {
      group.scale <- ggplot2::scale_colour_discrete()
    }
  }
  if (!inherits(group.scale, "Scale")) {
    stop("`group.scale` must be a ggplot2 colour or fill scale.", call. = FALSE)
  }
  group.scale <- group.scale$clone()
  group.aes <- group.scale$aesthetics
  if (!inherits(group.scale, if (continuous) "ScaleContinuous" else "ScaleDiscrete") ||
      length(group.aes) != 1L ||
      !group.aes %in% c("colour", "fill")) {
    stop("Use a ", if (continuous) "continuous" else "discrete",
         " colour or fill scale for group.scale.", call. = FALSE)
  }
  if (!continuous && !original.ids && inherits(group.scale$labels, "waiver")) {
    group.scale$labels <- function(v) {
      ifelse(is.na(v), "NA", paste(match(v, keys), v, sep = "  "))
    }
  }
  if (inherits(group.scale$name, "waiver")) {
    group.scale$name <- if (continuous) group.by else NULL
  }
  group.guide <- if (continuous) ggplot2::guide_colourbar(
    direction = "horizontal", order = 1,
    theme = ggplot2::theme(legend.title.position = "top",
                           legend.key.width = grid::unit(38, "mm"),
                           legend.key.height = grid::unit(1.8, "mm"),
                           legend.frame = ggplot2::element_blank())
  ) else ggplot2::guide_legend(
    ncol = min(4L, length(keys)), direction = "horizontal",
    override.aes = list(size = 2.0),
    theme = ggplot2::theme(legend.byrow = TRUE)
  )
  bar.guide <- ggplot2::guide_colourbar(
    position = if (continuous) "bottom" else "inside",
    direction = "horizontal", title = "Expression level", order = 2,
    theme = ggplot2::theme(legend.title.position = "top",
                           legend.key.width = grid::unit(83, "mm"),
                           legend.key.height = grid::unit(1.8, "mm"),
                           legend.frame = ggplot2::element_blank())
  )
  if (group.aes == "colour") {
    ref.layer <- ggplot2::geom_point(data = ref, ggplot2::aes(.data$x, .data$y, colour = .data$group),
                                     size = pt.size, stroke = 0)
    gene.layer <- ggplot2::geom_point(data = dots, ggplot2::aes(.data$x, .data$y, fill = .data$value),
                                      size = pt.size, shape = 21, colour = "transparent", stroke = 0)
    gradient <- ggplot2::scale_fill_gradientn(colours = feature.cols, guide = bar.guide)
    guides <- ggplot2::guides(colour = group.guide)
  } else {
    ref.layer <- ggplot2::geom_point(data = ref, ggplot2::aes(.data$x, .data$y, fill = .data$group),
                                     size = pt.size, shape = 21, colour = "transparent", stroke = 0)
    gene.layer <- ggplot2::geom_point(data = dots, ggplot2::aes(.data$x, .data$y, colour = .data$value),
                                      size = pt.size, stroke = 0)
    gradient <- ggplot2::scale_colour_gradientn(colours = feature.cols, guide = bar.guide)
    guides <- ggplot2::guides(fill = group.guide)
  }

  sep <- max(ref$x) + 0.13
  p <- ggplot2::ggplot() +
    ggplot2::annotate("segment", x = sep, xend = sep, y = bottom,
                      yend = title.y, colour = "#D9DFE5", linewidth = 0.28) +
    ggplot2::geom_polygon(data = frames, ggplot2::aes(.data$x, .data$y, group = .data$panel),
                          fill = NA, colour = "#D4DDE6", linewidth = 0.27) +
    ref.layer + gene.layer +
    ggplot2::geom_label(data = centres, ggplot2::aes(.data$x, .data$y, label = .data$id),
                        size = 2.7, fontface = "bold", colour = text.color,
                        label.padding = grid::unit(0.055, "lines"),
                        label.r = grid::unit(0.12, "lines"),
                        fill = "#FFFFFFBD", label.size = 0) +
    ggplot2::geom_text(data = titles, ggplot2::aes(.data$x, .data$y, label = .data$label),
                       hjust = 0, size = 3.6, fontface = "bold", colour = text.color) +
    # Equal-length perpendicular arms; coord_fixed preserves their screen lengths.
    ggplot2::annotate("segment", x = axis.x, xend = c(axis.x, axis.x + axis.length),
                      y = axis.y, yend = c(axis.y + axis.length, axis.y),
                      colour = "black", linewidth = 0.45, lineend = "butt") +
    ggplot2::annotate("text", x = axis.x + axis.length / 2, y = axis.y - 0.018,
                      label = embedding.label, vjust = 1, size = 3.2, colour = text.color) +
    group.scale + gradient + guides +
    ggplot2::coord_fixed(xlim = c(left, right),
                         ylim = c(bottom - 0.30, title.y + 0.06),
                         expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = text.color),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(7, 8, 5, 8),
      legend.position = "bottom",
      legend.justification.bottom = if (continuous) c(0, 0) else c(0.5, 0),
      legend.box = if (continuous) "horizontal" else "vertical",
      legend.position.inside = c(((gene.left + gene.right) / 2 - left) / (right - left), 0.003),
      legend.justification.inside = c(0.5, 0),
      legend.text = ggplot2::element_text(size = 7.8, colour = text.color),
      legend.title = ggplot2::element_text(size = 8.3, colour = text.color),
      legend.ticks = ggplot2::element_line(colour = text.color, linewidth = 0.3,
                                           lineend = "butt"),
      legend.ticks.length = grid::unit(c(-0.8, 0), "mm"),
      legend.key.width = grid::unit(3.3, "mm"),
      legend.key.height = grid::unit(3.8, "mm"),
      legend.spacing.x = grid::unit(if (continuous) 34 else 1.6, "mm"),
      legend.key.spacing.x = grid::unit(4, "mm"),
      legend.margin = ggplot2::margin(1, 0, 1, 0),
      legend.box.spacing = grid::unit(2, "mm")
    )
  p
}
