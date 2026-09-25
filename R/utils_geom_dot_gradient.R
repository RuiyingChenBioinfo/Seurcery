#' Add dot summaries below violin plots
#'
#' Creates a lightweight layer-like object that can be added to a ggplot object
#' with `+`. The corresponding S3 method `ggplot_add.geom_dot_gradient()`
#' extracts plotting data from an existing violin plot, calculates expressed
#' percentage and average expression for each identity and split group within
#' each facet, and draws summary dots beside the low end of the expression
#' axis, with optional significance labels and
#' an optional side guide panel.
#'
#' This function supports both a single ggplot violin plot and a patchwork
#' object such as that returned by `Seurat::VlnPlot()` when multiple features
#' are combined.
#'
#' When `split.by` is not used in `Seurat::VlnPlot()`, Seurat often maps fill to
#' the same grouping variable as x. In that case, this function automatically
#' collapses the comparison structure to one dot per x group, but the dot color
#' still follows x groups so each identity can have a different color.
#' Stacked Seurat plots retain feature colors and normalize average-expression
#' alpha separately within each feature panel. Dot sizes use one common
#' percentage scale across all panels of a faceted plot. Horizontal violins,
#' transformed scales, and existing coordinates are preserved. The function
#' can be combined with [geom_vln_gradient()] in either order.
#'
#' @param expr_cutoff Expression cutoff used to define whether a value is
#'   expressed.
#' @param group_cols Optional vector of colors. When `split.by` is used, its
#'   length must equal the number of split groups. When `split.by` is not used,
#'   its length must equal the number of x groups.
#' @param alpha_max Maximum alpha value.
#' @param alpha_min Minimum alpha value.
#' @param low_quantile Lower anchor point of the alpha interpolation, in
#'   `[0, 1]`.
#' @param mid_quantile Middle anchor point of the alpha interpolation, in
#'   `[0, 1]`.
#' @param high_quantile Upper anchor point of the alpha interpolation, in
#'   `[0, 1]`.
#' @param size_range Numeric vector of length 2 giving the dot size range for
#'   non-zero expressed percentages.
#' @param size_zero Dot size used when expressed percentage is zero.
#' @param size_power Power exponent used for nonlinear size scaling.
#' @param layer Optional index of the violin layer in the ggplot object. If
#'   `NULL`, the function tries to detect the first violin-like layer
#'   automatically.
#' @param dot_shape Shape passed to `ggplot2::geom_point()`.
#' @param dot_stroke Stroke width for the summary dots.
#' @param dodge_width Horizontal spread used to separate split groups.
#' @param star_size Text size for significance stars.
#' @param sig_test_all Optional test against all remaining observations. Must be
#'   `NULL`, `"wilcox"`, or `"ttest"`.
#' @param sig_adjust Adjustment method passed to `stats::p.adjust()`.
#' @param sig_cutoffs Named numeric vector defining p value cutoffs for
#'   significance stars.
#' @param guide Logical. Whether to draw the side guide panel.
#' @param guide_width Relative width of the guide panel in the final patchwork
#'   layout.
#' @param guide_size_breaks Breaks, in percent, to show in the expressed
#'   percentage size guide.
#' @param guide_title_size Text size for guide titles.
#' @param guide_label_size Text size for guide labels.
#' @param guide_bar_n Number of tiles used to draw the average expression guide
#'   bar.
#' @param guide_bar_low Low color for the average expression guide bar.
#' @param guide_bar_high High color for the average expression guide bar.
#' @param y_axis_min Requested lower bound of the expression axis (including
#'   horizontal violins). Default is `-1`. If this bound falls inside the data
#'   range or outside the scale's domain, space for the dots is added
#'   automatically below the displayed expression range. NULL also uses
#'   automatic placement.
#' @param min_avg_pct_diff Minimum absolute difference in expressed proportion
#'   required for significance stars to be shown. Uses proportion scale in
#'   `[0, 1]`. Default `0.2` means at least 20 percentage points.
#'
#' @return An object of class `"geom_dot_gradient"` that can be added to a
#'   ggplot object.
#'
#' @export
geom_dot_gradient <- function(
  expr_cutoff = 0,
  group_cols = NULL,
  alpha_max = 1,
  alpha_min = 0,
  low_quantile = 0,
  mid_quantile = 0.5,
  high_quantile = 1,
  size_range = c(0.03, 8),
  size_zero = 0.001,
  size_power = 2.2,
  layer = NULL,
  dot_shape = 21,
  dot_stroke = 0.15,
  dodge_width = 0.90,
  star_size = 4.5,
  sig_test_all = NULL,
  sig_adjust = "BH",
  sig_cutoffs = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
  guide = TRUE,
  guide_width = 0.42,
  guide_size_breaks = c(0, 20, 40, 60, 80, 100),
  guide_title_size = 12,
  guide_label_size = 10,
  guide_bar_n = 200,
  guide_bar_low = "grey97",
  guide_bar_high = "grey15",
  y_axis_min = -1,
  min_avg_pct_diff = 0.2
) {
  structure(
    list(
      expr_cutoff = expr_cutoff,
      group_cols = group_cols,
      alpha_max = alpha_max,
      alpha_min = alpha_min,
      low_quantile = low_quantile,
      mid_quantile = mid_quantile,
      high_quantile = high_quantile,
      size_range = size_range,
      size_zero = size_zero,
      size_power = size_power,
      layer = layer,
      dot_shape = dot_shape,
      dot_stroke = dot_stroke,
      dodge_width = dodge_width,
      star_size = star_size,
      sig_test_all = sig_test_all,
      sig_adjust = sig_adjust,
      sig_cutoffs = sig_cutoffs,
      guide = guide,
      guide_width = guide_width,
      guide_size_breaks = guide_size_breaks,
      guide_title_size = guide_title_size,
      guide_label_size = guide_label_size,
      guide_bar_n = guide_bar_n,
      guide_bar_low = guide_bar_low,
      guide_bar_high = guide_bar_high,
      y_axis_min = y_axis_min,
      min_avg_pct_diff = min_avg_pct_diff
    ),
    class = "geom_dot_gradient"
  )
}

#' @importFrom ggplot2 ggplot_add
#' @importFrom rlang .data
#' @method ggplot_add geom_dot_gradient
#' @export
ggplot_add.geom_dot_gradient <- function(object, plot, object_name) {
  .geom_dot_gradient_validate_args(object)
  .seurcery_apply_gradient(plot, object, .geom_dot_gradient_add_single_plot)
}

# internal helpers ---------------------------------------------------------

.geom_dot_gradient_interp_alpha_piecewise <- function(
  y,
  q_low_y,
  q_mid_y,
  q_high_y,
  alpha_low_end,
  alpha_mid,
  alpha_high_end
) {
  eps <- 1e-12

  if (y <= q_low_y) {
    return(alpha_low_end)
  }
  if (y >= q_high_y) {
    return(alpha_high_end)
  }
  if (abs(y - q_mid_y) <= eps) {
    return(alpha_mid)
  }
  if (y < q_mid_y) {
    if ((q_mid_y - q_low_y) <= eps) {
      return(alpha_mid)
    }
    return(
      alpha_low_end +
        (y - q_low_y) / (q_mid_y - q_low_y) * (alpha_mid - alpha_low_end)
    )
  }
  if ((q_high_y - q_mid_y) <= eps) {
    return(alpha_high_end)
  }

  alpha_mid +
    (y - q_mid_y) / (q_high_y - q_mid_y) * (alpha_high_end - alpha_mid)
}

.geom_dot_gradient_p_to_star <- function(p, cutoffs) {
  if (!is.finite(p) || is.na(p)) {
    return("")
  }

  labs <- names(cutoffs)
  vals <- unname(cutoffs)
  ord <- order(vals)
  vals <- vals[ord]
  labs <- labs[ord]

  hit <- which(p <= vals)
  if (length(hit) == 0) {
    return("")
  }

  labs[min(hit)]
}

.geom_dot_gradient_choose_display_max_break <- function(obs_max_pct, breaks_pct) {
  breaks_pct <- sort(unique(as.numeric(breaks_pct)))
  pos_breaks <- breaks_pct[breaks_pct > 0]

  if (length(pos_breaks) == 0) {
    return(100)
  }

  reached <- pos_breaks[pos_breaks >= obs_max_pct - 1e-12]
  if (length(reached) == 0) {
    return(max(100, obs_max_pct))
  }
  min(reached)
}

.geom_dot_gradient_size_from_pct <- function(
  pct_prop,
  obs_max_prop,
  size_range,
  size_zero,
  size_power
) {
  if (!is.finite(pct_prop) || is.na(pct_prop) || pct_prop <= 0) {
    return(size_zero)
  }
  if (!is.finite(obs_max_prop) || is.na(obs_max_prop) || obs_max_prop <= 0) {
    return(size_zero)
  }

  scaled_prop <- pct_prop / obs_max_prop
  scaled_prop <- max(min(scaled_prop, 1), 0)

  size_range[1] +
    (scaled_prop ^ size_power) * (size_range[2] - size_range[1])
}

.geom_dot_gradient_validate_args <- function(object) {
  .seurcery_validate_alpha(object)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  nonnegative <- c("size_zero", "dot_stroke", "dodge_width", "star_size",
                   "guide_title_size", "guide_label_size", "min_avg_pct_diff")
  for (nm in nonnegative) {
    if (!scalar(object[[nm]]) || object[[nm]] < 0) {
      stop("`", nm, "` must be a finite nonnegative number.", call. = FALSE)
    }
  }
  if (object$min_avg_pct_diff > 1) {
    stop("`min_avg_pct_diff` must be in [0, 1].", call. = FALSE)
  }
  for (nm in c("size_power", "guide_width")) {
    if (!scalar(object[[nm]]) || object[[nm]] <= 0) {
      stop("`", nm, "` must be a finite positive number.", call. = FALSE)
    }
  }
  if (!scalar(object$expr_cutoff)) {
    stop("`expr_cutoff` must be a finite number.", call. = FALSE)
  }
  if (!is.numeric(object$size_range) || length(object$size_range) != 2L ||
      any(!is.finite(object$size_range)) || object$size_range[1] < 0 ||
      object$size_range[2] <= object$size_range[1]) {
    stop("`size_range` must contain two finite increasing nonnegative numbers.", call. = FALSE)
  }
  if (!is.null(object$y_axis_min) && !scalar(object$y_axis_min)) {
    stop("`y_axis_min` must be NULL or a finite number.", call. = FALSE)
  }
  if (!is.logical(object$guide) || length(object$guide) != 1L || is.na(object$guide)) {
    stop("`guide` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!scalar(object$guide_bar_n) || object$guide_bar_n < 2 ||
      object$guide_bar_n != floor(object$guide_bar_n)) {
    stop("`guide_bar_n` must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(object$guide_size_breaks) || !length(object$guide_size_breaks) ||
      any(!is.finite(object$guide_size_breaks)) ||
      any(object$guide_size_breaks < 0 | object$guide_size_breaks > 100)) {
    stop("`guide_size_breaks` must contain finite percentages in [0, 100].", call. = FALSE)
  }
  if (!is.null(object$sig_test_all) &&
      (length(object$sig_test_all) != 1L ||
       !object$sig_test_all %in% c("wilcox", "ttest"))) {
    stop('`sig_test_all` must be NULL, "wilcox", or "ttest".', call. = FALSE)
  }
  if (length(object$sig_adjust) != 1L ||
      !object$sig_adjust %in% stats::p.adjust.methods) {
    stop("`sig_adjust` must be a method supported by p.adjust().", call. = FALSE)
  }
  if (!is.numeric(object$sig_cutoffs) || !length(object$sig_cutoffs) ||
      any(!is.finite(object$sig_cutoffs)) ||
      any(object$sig_cutoffs < 0 | object$sig_cutoffs > 1) ||
      is.null(names(object$sig_cutoffs)) || any(!nzchar(names(object$sig_cutoffs)))) {
    stop("`sig_cutoffs` must be a named numeric vector of cutoffs in [0, 1].", call. = FALSE)
  }
  invisible(TRUE)
}

# Evaluate the effective layer mappings, including .data[[...]], expressions,
# and layer-specific data. PANEL assignment comes from the trained facet.
.geom_dot_gradient_extract_plot_data <- function(plot_single, violin_layer, gb) {
  source <- plot_single$layers[[violin_layer]]
  raw_df <- source$data
  if (is.null(raw_df) || inherits(raw_df, "waiver")) raw_df <- plot_single$data
  if (is.function(raw_df)) raw_df <- raw_df(plot_single$data)
  if (!is.data.frame(raw_df)) {
    stop("The violin layer must use data-frame data.", call. = FALSE)
  }
  mapping <- if (isTRUE(source$inherit.aes)) plot_single$mapping else ggplot2::aes()
  mapping[names(source$mapping)] <- source$mapping
  evaluate <- function(nm, default = NULL) {
    if (is.null(mapping[[nm]])) return(default)
    value <- rlang::eval_tidy(mapping[[nm]], data = raw_df)
    if (length(value) == 1L) value <- rep(value, nrow(raw_df))
    if (length(value) != nrow(raw_df)) {
      stop("The `", nm, "` mapping does not match the violin data.", call. = FALSE)
    }
    value
  }
  xx <- evaluate("x")
  yy <- evaluate("y")
  if (is.null(xx) || is.null(yy)) {
    stop("Cannot evaluate the violin's x and y mappings.", call. = FALSE)
  }
  vdat <- gb$data[[violin_layer]]
  horizontal <- if ("flipped_aes" %in% names(vdat) && nrow(vdat)) {
    isTRUE(vdat$flipped_aes[1])
  } else {
    is.numeric(xx) && !is.numeric(yy)
  }
  group <- if (horizontal) yy else xx
  expr <- if (horizontal) xx else yy
  if (!is.numeric(expr)) {
    stop("The expression mapping must evaluate to numeric values.", call. = FALSE)
  }
  fill <- evaluate("fill", rep("All", nrow(raw_df)))
  raw_df$.x <- group
  raw_df$.y <- expr
  raw_df$.fill <- fill
  raw_df$.source_fill <- fill
  raw_df <- gb$layout$facet$map_data(
    raw_df, gb$layout$layout, gb$layout$facet_params
  )
  raw_df <- raw_df[is.finite(raw_df$.y) & !is.na(raw_df$.x) &
                     !is.na(raw_df$.fill), , drop = FALSE]
  if (!nrow(raw_df)) stop("No finite observations found for dot summaries.", call. = FALSE)
  facet_vars <- setdiff(names(gb$layout$layout),
                       c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y", "COORD"))
  facet_specs <- c(plot_single$facet$params$rows, plot_single$facet$params$cols,
                   plot_single$facet$params$facets)
  input_vars <- unique(unlist(lapply(facet_specs, function(q) {
    all.vars(rlang::get_expr(q))
  })))
  list(dat = raw_df, expression_axis = if (horizontal) "x" else "y",
       facet_vars = unique(c(facet_vars, intersect(input_vars, names(raw_df)))),
       has_fill_mapping = !is.null(mapping$fill),
       fixed_fill = source$aes_params$fill, fill_identity = inherits(fill, "AsIs"),
       fill_levels = if (is.factor(fill)) levels(fill) else unique(fill))
}

.geom_dot_gradient_levels <- function(x) {
  values <- if (is.factor(x)) levels(x) else unique(as.character(x))
  values[values %in% as.character(x)]
}

.geom_dot_gradient_add_facets <- function(data, panel, panel_data, extracted, gb) {
  if (is.null(data) || !nrow(data)) return(data)
  layout_row <- gb$layout$layout[as.character(gb$layout$layout$PANEL) == panel, , drop = FALSE]
  for (nm in extracted$facet_vars) {
    value <- if (nm %in% names(layout_row)) layout_row[[nm]][1] else panel_data[[nm]][1]
    data[[nm]] <- rep(value, nrow(data))
  }
  data
}

# Positions are computed in the trained expression scale then returned to data
# space. This avoids resetting log/reverse scales or coord_flip().
.geom_dot_gradient_positions <- function(panel_data, scale, object) {
  trans <- if (is.function(scale$get_transformation)) scale$get_transformation() else scale$trans
  values <- suppressWarnings(trans$transform(panel_data$.y))
  values <- values[is.finite(values)]
  if (!length(values)) stop("No observations lie in the expression scale's domain.", call. = FALSE)
  rng <- range(values)
  identity <- identical(trans$name, "identity")
  baseline <- if (identity) min(0, rng[1]) else rng[1]
  requested <- if (is.null(object$y_axis_min)) NA_real_ else
    suppressWarnings(trans$transform(object$y_axis_min))
  padding <- if (identity) max(1, diff(rng) * 0.1) else max(0.1, diff(rng) * 0.15)
  lower <- if (length(requested) == 1L && is.finite(requested) && requested < baseline) {
    requested
  } else {
    baseline - padding
  }
  list(dot = trans$inverse((baseline + lower) / 2),
       anchor = trans$inverse(lower), lower_transformed = lower)
}

# Preserve an explicitly limited expression scale while making room for dots.
# ggplot2 stores numeric continuous limits in transformed coordinates.
.geom_dot_gradient_extend_limits <- function(plot, axis, lower) {
  scale <- plot$scales$get_scales(axis)
  if (is.null(scale) || is.null(scale$limits)) return(plot)
  scale <- scale$clone()
  if (is.numeric(scale$limits)) {
    scale$limits[1] <- min(c(scale$limits[1], lower), na.rm = TRUE)
  } else if (is.function(scale$limits)) {
    old_limits <- scale$limits
    trans <- if (is.function(scale$get_transformation)) scale$get_transformation() else scale$trans
    scale$limits <- local({
      previous <- old_limits
      transform <- trans
      bound <- lower
      function(x) {
        lim <- transform$transform(previous(x))
        lim[1] <- min(c(lim[1], bound), na.rm = TRUE)
        transform$inverse(lim)
      }
    })
  }
  suppressMessages(plot + scale)
}

.geom_dot_gradient_extend_coord <- function(plot, axis, anchors) {
  limits <- plot$coordinates$limits
  bound <- limits[[axis]]
  if (is.numeric(bound) && length(bound) == 2L) {
    if (is.finite(bound[1])) bound[1] <- min(c(bound[1], anchors), na.rm = TRUE)
    if (is.finite(bound[2])) bound[2] <- max(c(bound[2], anchors), na.rm = TRUE)
    limits[[axis]] <- bound
    original_coord <- plot$coordinates
    plot$coordinates <- ggplot2::ggproto(NULL, original_coord, limits = limits)
  }
  plot
}

.geom_dot_gradient_build_summary_df <- function(dat, x_levels, fill_levels, object, collapsed_fill) {
  sum_df <- stats::aggregate(
    dat$.y,
    by = list(.x = dat$.x, .fill = dat$.fill),
    FUN = function(z) c(
      avg = mean(z, na.rm = TRUE),
      pct = mean(z > object$expr_cutoff, na.rm = TRUE) * 100,
      n = length(z)
    )
  )

  sum_df <- do.call(data.frame, sum_df)
  colnames(sum_df)[3:5] <- c("avg_expr", "pct_expr", "n")

  sum_df$.x <- factor(sum_df$.x, levels = x_levels)
  sum_df$.fill <- factor(sum_df$.fill, levels = fill_levels)
  sum_df <- sum_df[order(sum_df$.x, sum_df$.fill), , drop = FALSE]

  avg_rng <- range(sum_df$avg_expr, na.rm = TRUE)
  avg_den <- avg_rng[2] - avg_rng[1]
  if (!is.finite(avg_den) || avg_den <= 0) {
    avg_den <- 1
  }

  q_low_y <- avg_rng[1] + object$low_quantile * avg_den
  q_mid_y <- avg_rng[1] + object$mid_quantile * avg_den
  q_high_y <- avg_rng[1] + object$high_quantile * avg_den

  alpha_low_end <- object$alpha_min
  alpha_mid <- (object$alpha_min + object$alpha_max) / 2
  alpha_high_end <- object$alpha_max

  sum_df$alpha_val <- vapply(
    sum_df$avg_expr,
    function(y) {
      .geom_dot_gradient_interp_alpha_piecewise(
        y = y,
        q_low_y = q_low_y,
        q_mid_y = q_mid_y,
        q_high_y = q_high_y,
        alpha_low_end = alpha_low_end,
        alpha_mid = alpha_mid,
        alpha_high_end = alpha_high_end
      )
    },
    numeric(1)
  )

  obs_max_pct <- max(sum_df$pct_expr, na.rm = TRUE)
  display_max_pct <- .geom_dot_gradient_choose_display_max_break(
    obs_max_pct,
    object$guide_size_breaks
  )
  display_max_prop <- display_max_pct / 100

  sum_df$size_val <- vapply(
    sum_df$pct_expr / 100,
    function(p) {
      .geom_dot_gradient_size_from_pct(
        pct_prop = p,
        obs_max_prop = display_max_prop,
        size_range = object$size_range,
        size_zero = object$size_zero,
        size_power = object$size_power
      )
    },
    numeric(1)
  )

  if (isTRUE(collapsed_fill)) {
    sum_df$.disp_group <- as.character(sum_df$.x)
  } else {
    sum_df$.disp_group <- as.character(sum_df$.fill)
  }

  list(
    sum_df = sum_df,
    obs_max_pct = obs_max_pct,
    display_max_pct = display_max_pct,
    display_max_prop = display_max_prop
  )
}


.geom_dot_gradient_build_star_df <- function(
    dat,
    x_levels,
    fill_levels,
    collapsed_fill,
    star_base_y,
    object
) {
  if (is.null(object$sig_test_all)) {
    return(NULL)
  }

  if (isTRUE(collapsed_fill) || length(fill_levels) <= 1) {
    p_tab <- list()

    for (xx in x_levels) {
      g1 <- dat$.y[dat$.x == xx]
      g2 <- dat$.y[dat$.x != xx]

      if (length(g1) == 0 || length(g2) == 0) {
        next
      }

      pct1 <- mean(g1 > object$expr_cutoff, na.rm = TRUE)
      pct2 <- mean(g2 > object$expr_cutoff, na.rm = TRUE)
      pct_diff <- abs(pct1 - pct2)

      p_val <- tryCatch(
        {
          if (object$sig_test_all == "wilcox") {
            stats::wilcox.test(g1, g2, alternative = "greater")$p.value
          } else {
            stats::t.test(g1, g2, alternative = "greater")$p.value
          }
        },
        error = function(e) NA_real_
      )

      p_tab[[length(p_tab) + 1]] <- data.frame(
        .x = xx,
        .fill = "All",
        p = p_val,
        pct_diff = pct_diff,
        stringsAsFactors = FALSE
      )
    }

    if (length(p_tab) == 0) {
      return(NULL)
    }

    p_df <- do.call(rbind, p_tab)
    p_df$p_adj <- stats::p.adjust(p_df$p, method = object$sig_adjust)
    p_df$label <- vapply(
      p_df$p_adj,
      .geom_dot_gradient_p_to_star,
      character(1),
      cutoffs = object$sig_cutoffs
    )

    p_df <- p_df[
      nzchar(p_df$label) &
        is.finite(p_df$pct_diff) &
        !is.na(p_df$pct_diff) &
        p_df$pct_diff >= object$min_avg_pct_diff,
      ,
      drop = FALSE
    ]

    if (nrow(p_df) == 0) {
      return(NULL)
    }

    p_df$.x <- factor(p_df$.x, levels = x_levels)
    p_df$.fill <- factor(p_df$.fill, levels = fill_levels)
    p_df$y_plot <- star_base_y
    return(p_df)
  }

  p_tab <- list()

  for (xx in x_levels) {
    for (ff in fill_levels) {
      idx1 <- dat$.x == xx & dat$.fill == ff
      idx2 <- !idx1

      g1 <- dat$.y[idx1]
      g2 <- dat$.y[idx2]

      if (length(g1) == 0 || length(g2) == 0) {
        next
      }

      pct1 <- mean(g1 > object$expr_cutoff, na.rm = TRUE)
      pct2 <- mean(g2 > object$expr_cutoff, na.rm = TRUE)
      pct_diff <- abs(pct1 - pct2)

      p_val <- tryCatch(
        {
          if (object$sig_test_all == "wilcox") {
            stats::wilcox.test(g1, g2, alternative = "greater")$p.value
          } else {
            stats::t.test(g1, g2, alternative = "greater")$p.value
          }
        },
        error = function(e) NA_real_
      )

      p_tab[[length(p_tab) + 1]] <- data.frame(
        .x = xx,
        .fill = ff,
        p = p_val,
        pct_diff = pct_diff,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(p_tab) == 0) {
    return(NULL)
  }

  p_df <- do.call(rbind, p_tab)
  p_df$p_adj <- stats::p.adjust(p_df$p, method = object$sig_adjust)
  p_df$label <- vapply(
    p_df$p_adj,
    .geom_dot_gradient_p_to_star,
    character(1),
    cutoffs = object$sig_cutoffs
  )

  p_df <- p_df[
    nzchar(p_df$label) &
      is.finite(p_df$pct_diff) &
      !is.na(p_df$pct_diff) &
      p_df$pct_diff >= object$min_avg_pct_diff,
    ,
    drop = FALSE
  ]

  if (nrow(p_df) == 0) {
    return(NULL)
  }

  full_df <- expand.grid(
    .x = x_levels,
    .fill = fill_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  p_df <- merge(
    full_df,
    p_df,
    by = c(".x", ".fill"),
    all.x = TRUE,
    sort = FALSE
  )

  p_df$.x <- factor(p_df$.x, levels = x_levels)
  p_df$.fill <- factor(p_df$.fill, levels = fill_levels)
  p_df <- p_df[order(p_df$.x, p_df$.fill), , drop = FALSE]

  p_df$label[is.na(p_df$label)] <- ""
  p_df$y_plot <- star_base_y

  if (!any(nzchar(p_df$label))) {
    return(NULL)
  }

  return(p_df)
}

#' @keywords internal
#' @noRd
.geom_dot_gradient_add_single_plot <- function(plot_single, object) {
  violin_layer <- .seurcery_find_violin_layer(plot_single, object$layer)
  gb <- ggplot2::ggplot_build(plot_single)
  extracted <- .geom_dot_gradient_extract_plot_data(plot_single, violin_layer, gb)
  axis <- extracted$expression_axis
  group_axis <- if (axis == "y") "x" else "y"
  panel_ids <- unique(as.character(extracted$dat$PANEL))
  summaries <- stars <- anchors <- vector("list", length(panel_ids))
  lower_bounds <- numeric(length(panel_ids))
  global_groups <- .geom_dot_gradient_levels(extracted$dat$.x)
  global_fills <- .geom_dot_gradient_levels(extracted$dat$.fill)

  for (i in seq_along(panel_ids)) {
    panel <- panel_ids[i]
    dat <- extracted$dat[as.character(extracted$dat$PANEL) == panel, , drop = FALSE]
    scales <- gb$layout$get_scales(as.integer(panel))
    group_scale <- scales[[group_axis]]
    x_levels <- as.character(group_scale$get_limits())
    x_levels <- x_levels[x_levels %in% as.character(dat$.x)]
    if (!length(x_levels)) x_levels <- .geom_dot_gradient_levels(dat$.x)
    fill_levels <- .geom_dot_gradient_levels(dat$.fill)
    # A feature-colored stacked facet has one color per identity, but it is
    # not a split comparison. Keep its original fill for plotting below.
    collapsed <- all(vapply(split(as.character(dat$.fill), dat$.x, drop = TRUE),
                            function(z) length(unique(z)) <= 1L, logical(1)))
    dat$.x <- factor(dat$.x, levels = x_levels)
    dat$.fill <- factor(dat$.fill, levels = fill_levels)
    if (collapsed) {
      dat$.fill <- factor(rep("All", nrow(dat)), levels = "All")
      fill_levels <- "All"
    }
    result <- .geom_dot_gradient_build_summary_df(dat, x_levels, fill_levels,
                                                 object, collapsed)
    summary <- result$sum_df
    positions <- .geom_dot_gradient_positions(dat, scales[[axis]], object)
    lower_bounds[i] <- positions$lower_transformed
    summary$y_plot <- positions$dot
    summary$.source_fill <- if (collapsed) {
      as.character(dat$.source_fill[match(as.character(summary$.x), as.character(dat$.x))])
    } else {
      as.character(summary$.fill)
    }
    star <- .geom_dot_gradient_build_star_df(dat, x_levels, fill_levels,
                                            collapsed, positions$dot, object)
    center <- function(df) {
      value <- as.numeric(group_scale$map(as.character(df$.x)))
      if (!collapsed && length(fill_levels) > 1L) {
        value <- value + object$dodge_width *
          ((match(as.character(df$.fill), fill_levels) - 0.5) / length(fill_levels) - 0.5)
      }
      value
    }
    summary$.group_position <- center(summary)
    if (!is.null(star)) star$.group_position <- center(star)
    anchor <- data.frame(.group_position = as.numeric(group_scale$map(x_levels)),
                         y_plot = positions$anchor)

    if (!is.null(object$group_cols)) {
      keys <- if (collapsed) global_groups else global_fills
      cols <- object$group_cols
      if (is.null(names(cols))) {
        if (length(cols) != length(keys)) {
          stop("`group_cols` length must equal the number of displayed groups.", call. = FALSE)
        }
        names(cols) <- keys
      }
      if (!all(as.character(summary$.disp_group) %in% names(cols))) {
        stop("Named `group_cols` must contain all displayed group names.", call. = FALSE)
      }
      summary$.fixed_fill <- unname(cols[as.character(summary$.disp_group)])
    } else if (extracted$fill_identity) {
      summary$.fixed_fill <- summary$.source_fill
    } else if (!extracted$has_fill_mapping) {
      fill <- extracted$fixed_fill
      if (is.null(fill)) fill <- gb$data[[violin_layer]]$fill[1]
      if (!length(fill)) fill <- "grey70"
      summary$.fixed_fill <- rep(fill, nrow(summary))
    }
    summaries[[i]] <- .geom_dot_gradient_add_facets(summary, panel, dat, extracted, gb)
    stars[[i]] <- .geom_dot_gradient_add_facets(star, panel, dat, extracted, gb)
    anchors[[i]] <- .geom_dot_gradient_add_facets(anchor, panel, dat, extracted, gb)
  }
  sum_df <- do.call(rbind, summaries)
  star_df <- do.call(rbind, stars)
  anchor_df <- do.call(rbind, anchors)
  display_max_pct <- .geom_dot_gradient_choose_display_max_break(
    max(sum_df$pct_expr), object$guide_size_breaks
  )
  display_max_prop <- display_max_pct / 100
  sum_df$size_val <- vapply(sum_df$pct_expr / 100, .geom_dot_gradient_size_from_pct,
                          numeric(1), obs_max_prop = display_max_prop,
                          size_range = object$size_range, size_zero = object$size_zero,
                          size_power = object$size_power)
  if (extracted$has_fill_mapping && is.null(object$group_cols)) {
    sum_df$.source_fill <- factor(sum_df$.source_fill, levels = extracted$fill_levels)
  }
  position_mapping <- if (axis == "y") {
    ggplot2::aes(x = .data$.group_position, y = .data$y_plot)
  } else {
    ggplot2::aes(x = .data$y_plot, y = .data$.group_position)
  }
  dot_mapping <- position_mapping
  fixed_fill <- ".fixed_fill" %in% names(sum_df)
  if (!fixed_fill) dot_mapping$fill <- ggplot2::aes(fill = .data$.source_fill)$fill
  point_args <- list(data = sum_df, mapping = dot_mapping, inherit.aes = FALSE,
                     size = sum_df$size_val, alpha = sum_df$alpha_val,
                     shape = object$dot_shape, stroke = object$dot_stroke,
                     colour = "black", show.legend = FALSE)
  if (fixed_fill) point_args$fill <- sum_df$.fixed_fill
  main_plot <- plot_single +
    ggplot2::geom_blank(data = anchor_df, mapping = position_mapping, inherit.aes = FALSE) +
    do.call(ggplot2::geom_point, point_args)
  main_plot <- .geom_dot_gradient_extend_limits(main_plot, axis, min(lower_bounds))
  main_plot <- .geom_dot_gradient_extend_coord(main_plot, axis, anchor_df$y_plot)
  if (!is.null(star_df) && nrow(star_df)) {
    star_mapping <- position_mapping
    star_mapping$label <- ggplot2::aes(label = .data$label)$label
    main_plot <- main_plot + ggplot2::geom_text(
      data = star_df, mapping = star_mapping, inherit.aes = FALSE,
      colour = "black", size = object$star_size, vjust = 0.5,
      fontface = "plain", show.legend = FALSE
    )
  }
  if (!isTRUE(object$guide)) return(main_plot)

  size_breaks_pct <- sort(unique(as.numeric(object$guide_size_breaks)))
  size_breaks_pct <- size_breaks_pct[size_breaks_pct >= 0]
  if (length(size_breaks_pct) == 0) {
    size_breaks_pct <- c(0, 20, 40, 60, 80, 100)
  }

  size_breaks_pct <- size_breaks_pct[size_breaks_pct <= display_max_pct + 1e-12]
  if (!0 %in% size_breaks_pct) {
    size_breaks_pct <- c(0, size_breaks_pct)
  }
  if (!display_max_pct %in% size_breaks_pct) {
    size_breaks_pct <- c(size_breaks_pct, display_max_pct)
  }
  size_breaks_pct <- sort(unique(size_breaks_pct))

  guide_size_df <- data.frame(
    pct = size_breaks_pct,
    size_val = vapply(
      size_breaks_pct / 100,
      function(p) {
        .geom_dot_gradient_size_from_pct(
          pct_prop = p,
          obs_max_prop = display_max_prop,
          size_range = object$size_range,
          size_zero = object$size_zero,
          size_power = object$size_power
        )
      },
      numeric(1)
    )
  )
  guide_size_df$y <- rev(seq_len(nrow(guide_size_df)))
  guide_size_df$x <- 1

  bar_y <- seq(0, 1, length.out = object$guide_bar_n)
  guide_bar_df <- data.frame(
    x = 1,
    y = bar_y,
    fill_val = bar_y
  )

  guide_plot_size <- ggplot2::ggplot(
    guide_size_df,
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$size_val),
      shape = object$dot_shape,
      stroke = object$dot_stroke,
      fill = "grey70",
      colour = "black"
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::geom_text(
      data = guide_size_df,
      ggplot2::aes(
        x = 1.65,
        y = .data$y,
        label = paste0(.data$pct, "%")
      ),
      inherit.aes = FALSE,
      hjust = 0,
      size = object$guide_label_size / 3
    ) +
    ggplot2::annotate(
      "text",
      x = 1.2,
      y = max(guide_size_df$y) + 1.1,
      label = "Expressed\npercentage",
      size = object$guide_title_size / 3#,
      #fontface = "bold"
    ) +
    ggplot2::xlim(0.6, 3.6) +
    ggplot2::ylim(0.5, max(guide_size_df$y) + 1.8) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, 12, 5.5, 12)
    )

  guide_plot_bar <- ggplot2::ggplot(
    guide_bar_df,
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data$fill_val),
      height = 1 / object$guide_bar_n,
      width = 0.28
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c(object$guide_bar_low, object$guide_bar_high),
      limits = c(0, 1)
    ) +
    ggplot2::annotate(
      "text",
      x = 1.15,
      y = 1.2, #1.16
      label = "Average\nexpression",
      size = object$guide_title_size / 3#,
      #fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = 1.25, #1.45
      y = 0.1,
      label = "Low",
      hjust = 0,
      size = object$guide_label_size / 3
    ) +
    ggplot2::annotate(
      "text",
      x = 1.25, #1.45
      y = 0.9, #1
      label = "High",
      hjust = 0,
      size = object$guide_label_size / 3
    ) +
    ggplot2::xlim(0.7, 2.25) +
    ggplot2::ylim(-0.05, 1.28) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(5.5, 12, 5.5, 12)
    )

  guide_plot <- patchwork::wrap_plots(
    guide_plot_size,
    guide_plot_bar,
    ncol = 1,
    heights = c(0.62, 0.38)
  )

  patchwork::wrap_plots(
    main_plot,
    guide_plot,
    ncol = 2,
    widths = c(1, object$guide_width)
  )
}
