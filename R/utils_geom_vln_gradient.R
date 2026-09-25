#' Add gradient fill polygons to violin plots
#'
#' Creates a lightweight layer-like object that can be added to a ggplot object
#' with `+`. The object is handled by the S3 method
#' `ggplot_add.geom_vln_gradient()`, which reconstructs violin polygons from an
#' existing violin layer and fills them with an alpha gradient along the
#' expression axis, independently within each panel.
#'
#' This function supports both a single ggplot violin plot and a patchwork
#' object returned by `Seurat::VlnPlot()` when multiple features are combined.
#' Stacked Seurat plots retain their gene facets, manual colors, and orientation.
#' Gradients are rendered after scales are trained, so fill scales may be added
#' before or after this function. The quantile arguments are fractions of the
#' panel's expression range on the plotted scale, not empirical quantiles.
#' Both addition orders with [geom_dot_gradient()] are supported, including
#' nested patchwork plots containing the dot guides.
#'
#' @param bin Number of expression-axis bins used to reconstruct the gradient fill.
#'   Must be a finite integer greater than or equal to 2.
#' @param alpha_max Maximum alpha value.
#' @param alpha_min Minimum alpha value.
#' @param direction Gradient direction. Use `1` for low to high alpha from low
#'   expression to high expression, and `-1` for the reverse.
#' @param layer Optional index of the violin layer in the ggplot object. If
#'   `NULL`, the function tries to detect the first violin-like layer
#'   automatically.
#' @param outline Logical. Whether to redraw violin outlines.
#' @param outline_size Optional outline width. If `NULL`, uses the linewidth
#'   in the source violin layer when available.
#' @param min_width_frac Minimum local width threshold, expressed as a fraction
#'   of the maximum violin width, below which thin polygon slices are skipped.
#' @param high_quantile Upper anchor point of the alpha interpolation, in
#'   `[0, 1]`.
#' @param mid_quantile Middle anchor point of the alpha interpolation, in
#'   `[0, 1]`.
#' @param low_quantile Lower anchor point of the alpha interpolation, in
#'   `[0, 1]`.
#'
#' @return An object of class `"geom_vln_gradient"` that can be added to a
#'   ggplot object.
#' @export
geom_vln_gradient <- function(
  bin = 50,
  alpha_max = 1,
  alpha_min = 0,
  direction = 1,
  layer = NULL,
  outline = TRUE,
  outline_size = NULL,
  min_width_frac = 0.02,
  high_quantile = 1,
  mid_quantile = 0.5,
  low_quantile = 0
) {
  structure(
    list(
      bin = bin,
      alpha_max = alpha_max,
      alpha_min = alpha_min,
      direction = direction,
      layer = layer,
      outline = outline,
      outline_size = outline_size,
      min_width_frac = min_width_frac,
      high_quantile = high_quantile,
      mid_quantile = mid_quantile,
      low_quantile = low_quantile
    ),
    class = "geom_vln_gradient"
  )
}

#' Add a gradient violin layer to a ggplot object
#'
#' S3 method for adding a `"geom_vln_gradient"` object to a ggplot.
#'
#' @param object An object created by [geom_vln_gradient()].
#' @param plot A ggplot object or a patchwork object.
#' @param object_name Not used.
#'
#' @return A ggplot object or a patchwork object with reconstructed violin polygons.
#' @importFrom ggplot2 ggplot_add
#' @method ggplot_add geom_vln_gradient
#' @export
ggplot_add.geom_vln_gradient <- function(object, plot, object_name) {
  .seurcery_validate_alpha(object)
  if (!.seurcery_scalar(object$bin) || object$bin < 2 ||
      object$bin != floor(object$bin)) {
    stop("`bin` must be a finite integer >= 2.", call. = FALSE)
  }
  if (!.seurcery_scalar(object$direction) || !object$direction %in% c(-1, 1)) {
    stop("`direction` must be 1 or -1.", call. = FALSE)
  }
  if (!.seurcery_scalar(object$min_width_frac) || object$min_width_frac < 0 ||
      object$min_width_frac >= 1) {
    stop("`min_width_frac` must be in [0, 1).", call. = FALSE)
  }
  if (!is.logical(object$outline) || length(object$outline) != 1L ||
      is.na(object$outline)) stop("`outline` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(object$outline_size) &&
      (!.seurcery_scalar(object$outline_size) || object$outline_size < 0)) {
    stop("`outline_size` must be NULL or a non-negative finite scalar.", call. = FALSE)
  }
  .seurcery_apply_gradient(plot, object, .seurcery_add_violin)
}

.seurcery_add_violin <- function(plot, object) {
  index <- .seurcery_find_violin_layer(plot, object$layer)
  source_layer <- plot$layers[[index]]
  original <- source_layer$geom$.seurcery_source_geom
  if (is.null(original)) original <- source_layer$geom
  split_violin <- inherits(original, "GeomSplitViolin")
  # Derive a new layer: editing a shared ggproto would also alter the input plot.
  layer <- ggplot2::ggproto(NULL, source_layer)
  layer$geom <- ggplot2::ggproto(
    "GeomSeurceryViolinGradient", original,
    .seurcery_source_geom = original,
    parameters = function(self, extra = FALSE) original$parameters(extra),
    draw_panel = function(self, data, panel_params, coord, ...) {
      slices <- .seurcery_violin_slices(data, object, split_violin)
      fill_grob <- if (nrow(slices)) {
        ggplot2::GeomPolygon$draw_panel(slices, panel_params, coord)
      } else grid::nullGrob()
      if (!object$outline || !nrow(data)) return(fill_grob)
      data$fill <- NA
      data$alpha <- 1
      if (!is.null(object$outline_size)) {
        data$linewidth <- object$outline_size
        data$size <- object$outline_size
      }
      grid::grobTree(fill_grob, original$draw_panel(data, panel_params, coord, ...))
    }
  )
  plot$layers[[index]] <- layer
  plot
}

.seurcery_violin_slices <- function(data, object, split_violin = FALSE) {
  flipped <- "flipped_aes" %in% names(data) && isTRUE(data$flipped_aes[1L])
  data <- .seurcery_flip_data(data, flipped)
  needed <- c("x", "y", "xmin", "xmax", "violinwidth", "group")
  if (!all(needed %in% names(data))) {
    if (!nrow(data)) return(data.frame())
    stop("The violin layer lacks the density coordinates required for a gradient.", call. = FALSE)
  }
  data <- data[is.finite(data$y) & is.finite(data$x) &
                 is.finite(data$violinwidth), , drop = FALSE]
  if (!nrow(data)) return(data.frame())
  # draw_panel receives one facet at a time, so each gene has its own anchors.
  yrange <- range(data$y)
  if (diff(yrange) <= 0) return(data.frame())
  breaks <- seq(yrange[1L], yrange[2L], length.out = object$bin + 1L)
  mids <- (utils::head(breaks, -1L) + utils::tail(breaks, -1L)) / 2
  anchors <- yrange[1L] + c(object$low_quantile, object$mid_quantile,
                           object$high_quantile) * diff(yrange)
  ends <- c(object$alpha_min, object$alpha_max)
  if (object$direction == -1) ends <- rev(ends)
  alphas <- vapply(mids, .geom_dot_gradient_interp_alpha_piecewise, numeric(1),
                   q_low_y = anchors[1L], q_mid_y = anchors[2L],
                   q_high_y = anchors[3L], alpha_low_end = ends[1L],
                   alpha_mid = mean(ends), alpha_high_end = ends[2L])
  result <- list()
  id <- 0L
  for (g in split(data, data$group)) {
    left <- g$x - g$violinwidth * (g$x - g$xmin)
    right <- g$x + g$violinwidth * (g$xmax - g$x)
    curve <- stats::aggregate(cbind(left, right), list(y = g$y), mean)
    curve <- curve[order(curve$y), , drop = FALSE]
    if (nrow(curve) < 2L) next
    if (split_violin) {
      if (g$group[1L] %% 2L == 1L) curve$right <- g$x[1L]
      else curve$left <- g$x[1L]
    }
    cutoff <- max(curve$right - curve$left) * object$min_width_frac
    # Interpolate both boundaries once per violin instead of four times per bin.
    knots <- sort(unique(c(curve$y[1L], breaks[breaks > min(curve$y) &
                        breaks < max(curve$y)], curve$y[nrow(curve)])))
    xl <- stats::approx(curve$y, curve$left, knots, ties = mean)$y
    xr <- stats::approx(curve$y, curve$right, knots, ties = mean)$y
    for (i in seq_len(length(knots) - 1L)) {
      if (max(xr[i] - xl[i], xr[i + 1L] - xl[i + 1L]) <= cutoff) next
      id <- id + 1L
      bin_id <- findInterval(mean(knots[c(i, i + 1L)]), breaks,
                             all.inside = TRUE)
      result[[id]] <- data.frame(
        x = c(xl[i], xl[i + 1L], xr[i + 1L], xr[i]),
        y = knots[c(i, i + 1L, i + 1L, i)], group = id,
        fill = g$fill[1L], alpha = alphas[bin_id], colour = NA_character_,
        linewidth = 0, linetype = 1
      )
    }
  }
  if (!length(result)) return(data.frame())
  .seurcery_flip_data(do.call(rbind, result), flipped)
}
