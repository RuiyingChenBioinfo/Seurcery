#' Add dot summaries below violin plots
#'
#' Creates a lightweight layer-like object that can be added to a ggplot object
#' with `+`. The corresponding S3 method `ggplot_add.geom_dot_gradient()`
#' extracts plotting data from an existing violin plot, calculates expressed
#' percentage and average expression for each x and split group, and draws
#' summary dots below the violin panel, with optional significance labels and
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
#' @param y_axis_min Lower bound for the y axis. Default is `-1`.
#' @param min_avg_pct_diff Minimum absolute difference in expressed proportion
#'   required for significance stars to be shown. Uses proportion scale in
#'   `[0, 1]`. Default `0.1` means at least 10 percentage points.
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
  if (!inherits(plot, "ggplot") && !.geom_dot_gradient_is_patchwork_like(plot)) {
    stop("`geom_dot_gradient()` must be added to a ggplot or patchwork object.", call. = FALSE)
  }

  if (.geom_dot_gradient_is_patchwork_like(plot)) {
    return(.geom_dot_gradient_apply_to_patchwork(plot, object, .geom_dot_gradient_add_single_plot))
  }

  .geom_dot_gradient_add_single_plot(plot, object)
}

# internal helpers ---------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

.geom_dot_gradient_first_non_null <- function(...) {
  xs <- list(...)
  for (i in seq_along(xs)) {
    if (!is.null(xs[[i]]) && length(xs[[i]]) > 0) {
      return(xs[[i]])
    }
  }
  NULL
}

.geom_dot_gradient_aes_name <- function(aes_obj) {
  if (is.null(aes_obj) || length(aes_obj) == 0) {
    return(NULL)
  }

  expr <- tryCatch(
    {
      if (inherits(aes_obj, "quosure")) {
        rlang::get_expr(aes_obj)
      } else {
        aes_obj
      }
    },
    error = function(e) aes_obj
  )

  if (is.symbol(expr) || is.name(expr)) {
    return(as.character(expr))
  }

  if (is.character(expr) && length(expr) == 1) {
    return(expr)
  }

  NULL
}

.geom_dot_gradient_is_patchwork_like <- function(p) {
  inherits(p, "patchwork") && !is.null(p$patches)
}

.geom_dot_gradient_apply_to_patchwork <- function(pw, object, add_fun) {
  strip_patchwork <- function(x) {
    y <- x
    class(y) <- setdiff(class(y), "patchwork")
    y$patches <- NULL
    y
  }

  main_plot <- strip_patchwork(pw)
  all_plots <- c(list(main_plot), pw$patches$plots)
  all_plots <- lapply(all_plots, add_fun, object = object)

  lay <- pw$patches$layout
  ann <- pw$patches$annotation

  out <- patchwork::wrap_plots(
    all_plots,
    ncol = lay$ncol %||% NULL,
    nrow = lay$nrow %||% NULL,
    byrow = lay$byrow %||% NULL,
    guides = lay$guides %||% NULL
  )

  if (!is.null(lay$widths) || !is.null(lay$heights)) {
    out <- out + patchwork::plot_layout(
      widths = lay$widths %||% NULL,
      heights = lay$heights %||% NULL
    )
  }

  if (!is.null(ann)) {
    out <- out + patchwork::plot_annotation(
      title = ann$title %||% NULL,
      subtitle = ann$subtitle %||% NULL,
      caption = ann$caption %||% NULL,
      tag_levels = ann$tag_levels %||% NULL,
      tag_prefix = ann$tag_prefix %||% NULL,
      tag_suffix = ann$tag_suffix %||% NULL,
      tag_sep = ann$tag_sep %||% NULL,
      theme = ann$theme %||% NULL
    )
  }

  out
}

.geom_dot_gradient_find_violin_layer <- function(plot_obj, layer = NULL) {
  if (!is.null(layer)) {
    return(layer)
  }

  geom_names <- vapply(
    plot_obj$layers,
    function(z) paste(class(z$geom), collapse = "/"),
    character(1)
  )

  hit <- grep("SplitViolin|Violin", geom_names, ignore.case = TRUE)
  if (length(hit) == 0) {
    stop(
      "Cannot find a violin layer automatically. Please set `layer` manually.",
      call. = FALSE
    )
  }

  hit[1]
}

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

.geom_dot_gradient_default_seurat_cols <- function(n) {
  if (n <= 0) {
    return(character(0))
  }

  if (requireNamespace("Seurat", quietly = TRUE) &&
      "DiscretePalette" %in% getNamespaceExports("Seurat")) {
    cols <- tryCatch(
      Seurat::DiscretePalette(n = n),
      error = function(e) NULL
    )
    if (!is.null(cols) && length(cols) >= n) {
      return(unname(cols[seq_len(n)]))
    }
  }

  grDevices::hcl.colors(n = n, palette = "Dark 3")
}

#' @keywords internal
#' @noRd
.geom_dot_gradient_resolve_fill_cols <- function(
  group_cols,
  levels_to_color,
  gb = NULL,
  vdat = NULL,
  collapsed_fill = FALSE,
  x_levels = NULL
) {
  n_col <- length(levels_to_color)

  if (!is.null(group_cols)) {
    gc <- group_cols

    if (is.null(names(gc))) {
      if (length(gc) != n_col) {
        stop(
          "`group_cols` length must equal the number of displayed groups.",
          call. = FALSE
        )
      }
      names(gc) <- levels_to_color
    } else {
      if (!all(levels_to_color %in% names(gc))) {
        stop(
          "Named `group_cols` must contain all displayed group names.",
          call. = FALSE
        )
      }
      gc <- gc[levels_to_color]
      if (length(gc) != n_col) {
        stop(
          "`group_cols` length must equal the number of displayed groups.",
          call. = FALSE
        )
      }
    }

    gc <- unname(gc)
    names(gc) <- levels_to_color
    return(gc)
  }

  if (isTRUE(collapsed_fill) && !is.null(vdat) && !is.null(x_levels)) {
    if ("group" %in% colnames(vdat) && "fill" %in% colnames(vdat) && "x" %in% colnames(vdat)) {
      vg <- split(vdat, vdat$group)

      centers <- vapply(
        vg,
        function(df) stats::median(df$x, na.rm = TRUE),
        numeric(1)
      )
      fills <- vapply(
        vg,
        function(df) {
          vals <- unique(as.character(df$fill))
          vals <- vals[!is.na(vals) & nzchar(vals)]
          if (length(vals) == 0) NA_character_ else vals[1]
        },
        character(1)
      )

      ord <- order(centers)
      fills <- fills[ord]
      fills <- fills[seq_len(min(length(fills), length(x_levels)))]

      if (length(fills) == length(x_levels) && !any(is.na(fills))) {
        fills <- unname(fills)
        names(fills) <- x_levels
        return(fills)
      }
    }
  }

  if (!is.null(gb)) {
    fill_scale <- gb$plot$scales$get_scales("fill")
    if (!is.null(fill_scale)) {
      fill_cols <- tryCatch(
        fill_scale$map(levels_to_color),
        error = function(e) NULL
      )

      if (!is.null(fill_cols) &&
          length(fill_cols) == n_col &&
          !any(is.na(fill_cols))) {
        fill_cols <- unname(fill_cols)
        names(fill_cols) <- levels_to_color
        return(fill_cols)
      }

      fill_cols <- tryCatch(
        fill_scale$palette(n_col),
        error = function(e) NULL
      )

      if (!is.null(fill_cols) &&
          length(fill_cols) >= n_col &&
          !any(is.na(fill_cols[seq_len(n_col)]))) {
        fill_cols <- unname(fill_cols[seq_len(n_col)])
        names(fill_cols) <- levels_to_color
        return(fill_cols)
      }
    }
  }

  cols <- .geom_dot_gradient_default_seurat_cols(n_col)
  names(cols) <- levels_to_color
  cols
}

.geom_dot_gradient_choose_display_max_break <- function(obs_max_pct, breaks_pct) {
  breaks_pct <- sort(unique(as.numeric(breaks_pct)))
  pos_breaks <- breaks_pct[breaks_pct > 0]

  if (length(pos_breaks) == 0) {
    return(100)
  }

  reached <- pos_breaks[pos_breaks <= obs_max_pct + 1e-12]
  if (length(reached) == 0) {
    return(min(pos_breaks))
  }

  max(reached)
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
  q_vals <- c(object$low_quantile, object$mid_quantile, object$high_quantile)
  if (any(q_vals < 0) || any(q_vals > 1)) {
    stop(
      "`low_quantile`, `mid_quantile`, and `high_quantile` must be in [0, 1].",
      call. = FALSE
    )
  }

  if (!(object$low_quantile <= object$mid_quantile &&
        object$mid_quantile <= object$high_quantile)) {
    stop("Need low_quantile <= mid_quantile <= high_quantile.", call. = FALSE)
  }

  if (!is.numeric(object$size_range) ||
      length(object$size_range) != 2 ||
      object$size_range[1] < 0 ||
      object$size_range[2] <= object$size_range[1]) {
    stop(
      "`size_range` must be a numeric vector of length 2 with increasing values.",
      call. = FALSE
    )
  }

  if (!is.numeric(object$size_zero) ||
      length(object$size_zero) != 1 ||
      object$size_zero < 0) {
    stop("`size_zero` must be a non negative number.", call. = FALSE)
  }

  if (!is.numeric(object$size_power) ||
      length(object$size_power) != 1 ||
      object$size_power <= 0) {
    stop("`size_power` must be a positive number.", call. = FALSE)
  }

  if (!is.null(object$sig_test_all) &&
      !object$sig_test_all %in% c("wilcox", "ttest")) {
    stop(
      "`sig_test_all` must be NULL, \"wilcox\", or \"ttest\".",
      call. = FALSE
    )
  }

  if (!is.null(object$y_axis_min)) {
    if (!is.numeric(object$y_axis_min) ||
        length(object$y_axis_min) != 1 ||
        is.na(object$y_axis_min)) {
      stop("`y_axis_min` must be NULL or a single numeric value.", call. = FALSE)
    }
  }

  if (!is.numeric(object$min_avg_pct_diff) ||
      length(object$min_avg_pct_diff) != 1 ||
      is.na(object$min_avg_pct_diff) ||
      object$min_avg_pct_diff < 0) {
    stop("`min_avg_pct_diff` must be a single non negative numeric value.", call. = FALSE)
  }

  invisible(TRUE)
}

.geom_dot_gradient_extract_plot_data <- function(plot_single, violin_layer) {
  raw_df <- plot_single$data
  if (!is.data.frame(raw_df)) {
    stop(
      "Cannot extract plotting data from the current plot object.",
      call. = FALSE
    )
  }

  x_var <- .geom_dot_gradient_aes_name(.geom_dot_gradient_first_non_null(
    plot_single$mapping$x,
    plot_single$layers[[violin_layer]]$mapping$x
  ))
  y_var <- .geom_dot_gradient_aes_name(.geom_dot_gradient_first_non_null(
    plot_single$mapping$y,
    plot_single$layers[[violin_layer]]$mapping$y
  ))
  fill_var <- .geom_dot_gradient_aes_name(.geom_dot_gradient_first_non_null(
    plot_single$mapping$fill,
    plot_single$layers[[violin_layer]]$mapping$fill
  ))

  if (is.null(x_var) || is.null(y_var)) {
    stop("Cannot determine x and y variables from the current plot.", call. = FALSE)
  }

  if (!x_var %in% colnames(raw_df) || !y_var %in% colnames(raw_df)) {
    stop(
      "The x or y variables inferred from the plot are not found in plot$data.",
      call. = FALSE
    )
  }

  if (is.null(fill_var) || !fill_var %in% colnames(raw_df)) {
    fill_var <- NULL
  }

  keep_cols <- unique(c(x_var, y_var, fill_var))
  dat <- raw_df[, keep_cols, drop = FALSE]

  colnames(dat)[colnames(dat) == x_var] <- ".x"
  colnames(dat)[colnames(dat) == y_var] <- ".y"

  if (!is.null(fill_var) && fill_var %in% colnames(dat)) {
    colnames(dat)[colnames(dat) == fill_var] <- ".fill"
  }

  dat$.x <- as.character(dat$.x)
  collapsed_fill <- FALSE

  if (!(".fill" %in% colnames(dat))) {
    dat$.fill <- "All"
    fill_levels <- "All"
    collapsed_fill <- TRUE
  } else {
    dat$.fill <- as.character(dat$.fill)

    if (is.factor(raw_df[[fill_var]])) {
      fill_levels <- levels(raw_df[[fill_var]])
    } else {
      fill_levels <- unique(as.character(dat$.fill))
    }

    fill_levels <- fill_levels[fill_levels %in% unique(dat$.fill)]
    fill_levels <- fill_levels[!is.na(fill_levels) & nzchar(fill_levels)]

    if (length(fill_levels) == 0) {
      dat$.fill <- "All"
      fill_levels <- "All"
      collapsed_fill <- TRUE
    }
  }

  x_levels <- if (is.factor(raw_df[[x_var]])) {
    levels(raw_df[[x_var]])
  } else {
    unique(as.character(raw_df[[x_var]]))
  }
  x_levels <- x_levels[x_levels %in% unique(dat$.x)]

  dat <- dat[is.finite(dat$.y) & !is.na(dat$.x) & !is.na(dat$.fill), , drop = FALSE]
  if (nrow(dat) == 0) {
    stop("No valid data found for dot summary.", call. = FALSE)
  }

  if (!collapsed_fill) {
    same_as_x <- length(dat$.fill) == length(dat$.x) &&
      all(as.character(dat$.fill) == as.character(dat$.x))

    if (isTRUE(same_as_x)) {
      dat$.fill <- "All"
      fill_levels <- "All"
      collapsed_fill <- TRUE
    }
  }

  list(
    dat = dat,
    raw_df = raw_df,
    x_levels = x_levels,
    fill_levels = fill_levels,
    collapsed_fill = collapsed_fill
  )
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
            stats::t.test(g1, g2)$p.value
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
  
  p_df$.x <- factor(p_df$.x, levels = x_levels)
  p_df$.fill <- factor(p_df$.fill, levels = fill_levels)
  p_df$y_plot <- star_base_y
  p_df
}

#' @keywords internal
#' @noRd
.geom_dot_gradient_add_single_plot <- function(plot_single, object) {
  if (!inherits(plot_single, "ggplot")) {
    stop("`geom_dot_gradient()` must be added to a ggplot object.", call. = FALSE)
  }

  .geom_dot_gradient_validate_args(object)

  violin_layer <- .geom_dot_gradient_find_violin_layer(plot_single, object$layer)
  gb <- ggplot2::ggplot_build(plot_single)
  vdat <- gb$data[[violin_layer]]

  extracted <- .geom_dot_gradient_extract_plot_data(plot_single, violin_layer)
  dat <- extracted$dat
  x_levels <- extracted$x_levels
  fill_levels <- extracted$fill_levels
  collapsed_fill <- extracted$collapsed_fill

  summary_res <- .geom_dot_gradient_build_summary_df(
    dat = dat,
    x_levels = x_levels,
    fill_levels = fill_levels,
    object = object,
    collapsed_fill = collapsed_fill
  )
  sum_df <- summary_res$sum_df
  display_max_pct <- summary_res$display_max_pct
  display_max_prop <- summary_res$display_max_prop

  levels_to_color <- if (isTRUE(collapsed_fill)) x_levels else fill_levels
  fill_cols <- .geom_dot_gradient_resolve_fill_cols(
    group_cols = object$group_cols,
    levels_to_color = levels_to_color,
    gb = gb,
    vdat = vdat,
    collapsed_fill = collapsed_fill,
    x_levels = x_levels
  )

  sum_df$fill_col <- unname(fill_cols[sum_df$.disp_group])
  sum_df$fill_col_asis <- I(sum_df$fill_col)

  bottom_limit <- if (is.null(object$y_axis_min)) {
    -1
  } else {
    object$y_axis_min
  }

  if (bottom_limit <= -0.5) {
    dot_y <- -0.5
  } else {
    dot_y <- bottom_limit + 0.5 * (0 - bottom_limit)
  }

  star_base_y <- dot_y # + 0.06
  sum_df$y_plot <- dot_y

  star_df <- .geom_dot_gradient_build_star_df(
    dat = dat,
    x_levels = x_levels,
    fill_levels = fill_levels,
    collapsed_fill = collapsed_fill,
    star_base_y = star_base_y,
    object = object
  )

  anchor_df <- data.frame(
    .x = factor(x_levels, levels = x_levels),
    .y = rep(bottom_limit, length(x_levels))
  )

  if (isTRUE(collapsed_fill) || length(fill_levels) <= 1) {
    main_plot <- plot_single +
      ggplot2::geom_blank(
        data = anchor_df,
        mapping = ggplot2::aes(x = .data$.x, y = .data$.y),
        inherit.aes = FALSE
      ) +
      ggplot2::geom_point(
        data = sum_df,
        mapping = ggplot2::aes(
          x = .data$.x,
          y = .data$y_plot,
          size = .data$size_val,
          alpha = .data$alpha_val,
          fill = .data$fill_col_asis
        ),
        inherit.aes = FALSE,
        shape = object$dot_shape,
        stroke = object$dot_stroke,
        colour = "black",
        show.legend = FALSE
      ) +
      ggplot2::scale_size_identity() +
      ggplot2::scale_alpha_identity() +
      ggplot2::scale_y_continuous(limits = c(bottom_limit, NA)) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::guides(fill = "none", size = "none", alpha = "none") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(5.5, 10, 20, 5.5)
      )

    if (!is.null(star_df) && nrow(star_df) > 0) {
      main_plot <- main_plot +
        ggplot2::geom_text(
          data = star_df,
          mapping = ggplot2::aes(
            x = .data$.x,
            y = .data$y_plot,
            label = .data$label
          ),
          inherit.aes = FALSE,
          colour = "black",
          size = object$star_size,
          vjust = 0.5,
          fontface = "plain",
          show.legend = FALSE
        )
    }
  } else {
    dodge_pos <- ggplot2::position_dodge(width = object$dodge_width)

    main_plot <- plot_single +
      ggplot2::geom_blank(
        data = anchor_df,
        mapping = ggplot2::aes(x = .data$.x, y = .data$.y),
        inherit.aes = FALSE
      ) +
      ggplot2::geom_point(
        data = sum_df,
        mapping = ggplot2::aes(
          x = .data$.x,
          y = .data$y_plot,
          group = .data$.fill,
          size = .data$size_val,
          alpha = .data$alpha_val,
          fill = .data$fill_col_asis
        ),
        position = dodge_pos,
        inherit.aes = FALSE,
        shape = object$dot_shape,
        stroke = object$dot_stroke,
        colour = "black",
        show.legend = FALSE
      ) +
      ggplot2::scale_size_identity() +
      ggplot2::scale_alpha_identity() +
      ggplot2::scale_y_continuous(limits = c(bottom_limit, NA)) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::guides(fill = "none", size = "none", alpha = "none") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(5.5, 10, 20, 5.5)
      )

    if (!is.null(star_df) && nrow(star_df) > 0) {
      main_plot <- main_plot +
        ggplot2::geom_text(
          data = star_df,
          mapping = ggplot2::aes(
            x = .data$.x,
            y = .data$y_plot,
            label = .data$label,
            group = .data$.fill
          ),
          position = dodge_pos,
          inherit.aes = FALSE,
          colour = "black",
          size = object$star_size,
          vjust = 0.5,
          fontface = "plain",
          show.legend = FALSE
        )
    }
  }

  if (!isTRUE(object$guide)) {
    return(main_plot)
  }

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
      size = object$guide_title_size / 3,
      fontface = "bold"
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
      y = 1.09, #1.12
      label = "Average\nexpression",
      size = object$guide_title_size / 3,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = 1.45,
      y = 0,
      label = "Low",
      hjust = 0,
      size = object$guide_label_size / 3
    ) +
    ggplot2::annotate(
      "text",
      x = 1.45,
      y = 1,
      label = "High",
      hjust = 0,
      size = object$guide_label_size / 3
    ) +
    ggplot2::xlim(0.7, 2.25) +
    ggplot2::ylim(-0.05, 1.18) +
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