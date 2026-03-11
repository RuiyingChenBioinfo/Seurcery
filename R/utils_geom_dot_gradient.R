#' Add dot summaries with alpha gradients below violin plots
#'
#' Creates a lightweight layer-like object that can be added to a ggplot object
#' with `+`. The corresponding S3 method `ggplot_add.geom_dot_gradient()`
#' extracts plotting data from an existing violin plot, calculates expressed
#' percentage and average expression for each x and split group, and draws
#' summary dots with gradient alpha, optional significance labels, and an
#' optional side guide panel.
#'
#' @param expr_cutoff Expression cutoff used to define whether a value is
#'   expressed.
#' @param group_cols Optional vector of fill colors for split groups. Can be
#'   named or unnamed.
#' @param alpha_max Maximum alpha value.
#' @param alpha_min Minimum alpha value.
#' @param direction Gradient direction. Use `1` for low to high alpha from low
#'   average expression to high average expression, and `-1` for the reverse.
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
#' @param y_offset Vertical offset of summary dots below the violin panel,
#'   expressed as a fraction of y-range.
#' @param dodge_width Horizontal spread used to separate split groups.
#' @param star_offset Vertical offset for significance stars, expressed as a
#'   fraction of y-range.
#' @param star_size Text size for significance stars.
#' @param sig_test_all Optional test against all remaining observations. Must be
#'   `NULL`, `"wilcox"`, or `"ttest"`.
#' @param sig_adjust Adjustment method passed to `stats::p.adjust()`.
#' @param sig_cutoffs Named numeric vector defining p-value cutoffs for
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
#'
#' @return An object of class `"geom_dot_gradient"` that can be added to a
#'   ggplot object.
#' @export
geom_dot_gradient <- function(expr_cutoff = 0,
                              group_cols = NULL,
                              alpha_max = 1,
                              alpha_min = 0,
                              direction = 1,
                              low_quantile = 0,
                              mid_quantile = 0.5,
                              high_quantile = 1,
                              size_range = c(0.03, 8),
                              size_zero = 0.001,
                              size_power = 2.2,
                              layer = NULL,
                              dot_shape = 21,
                              dot_stroke = 0.15,
                              y_offset = 0.10,
                              dodge_width = 0.90,
                              star_offset = 0.012,
                              star_size = 4.5,
                              sig_test_all = NULL,
                              sig_adjust = "BH",
                              sig_cutoffs = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
                              guide = TRUE,
                              guide_width = 0.32,
                              guide_size_breaks = c(0, 20, 40, 60, 80, 100),
                              guide_title_size = 12,
                              guide_label_size = 10,
                              guide_bar_n = 200,
                              guide_bar_low = "grey97",
                              guide_bar_high = "grey15") {
  structure(
    list(
      expr_cutoff = expr_cutoff,
      group_cols = group_cols,
      alpha_max = alpha_max,
      alpha_min = alpha_min,
      direction = direction,
      low_quantile = low_quantile,
      mid_quantile = mid_quantile,
      high_quantile = high_quantile,
      size_range = size_range,
      size_zero = size_zero,
      size_power = size_power,
      layer = layer,
      dot_shape = dot_shape,
      dot_stroke = dot_stroke,
      y_offset = y_offset,
      dodge_width = dodge_width,
      star_offset = star_offset,
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
      guide_bar_high = guide_bar_high
    ),
    class = "geom_dot_gradient"
  )
}

#' Add a dot gradient layer to a ggplot object
#'
#' S3 method for adding a `"geom_dot_gradient"` object to a ggplot.
#'
#' @param object An object created by [geom_dot_gradient()].
#' @param plot A ggplot object.
#' @param object_name Not used.
#'
#' @return A ggplot object or a patchwork combined plot.
#' @export
ggplot_add.geom_dot_gradient <- function(object, plot, object_name) {
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0) y else x
  }

  aes_name <- function(aes_obj) {
    if (is.null(aes_obj)) {
      return(NULL)
    }

    expr <- tryCatch(
      rlang::get_expr(aes_obj),
      error = function(e) aes_obj
    )

    if (is.symbol(expr) || is.name(expr)) {
      return(as.character(expr))
    }
    if (is.character(expr) && length(expr) == 1) {
      return(expr)
    }

    txt <- paste(deparse(expr), collapse = "")
    txt <- gsub("^~", "", txt)
    txt
  }

  find_violin_layer <- function(plot, layer = NULL) {
    if (!is.null(layer)) {
      return(layer)
    }

    geom_names <- vapply(
      plot$layers,
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

  interp_alpha_piecewise <- function(y, q_low_y, q_mid_y, q_high_y,
                                     alpha_low_end, alpha_mid, alpha_high_end) {
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

  p_to_star <- function(p, cutoffs) {
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

  subgroup_offsets <- function(n, width) {
    if (n <= 1) {
      return(0)
    }
    width * (((seq_len(n) - 0.5) / n) - 0.5)
  }

  resolve_fill_cols <- function(group_cols, fill_levels, gb) {
    if (!is.null(group_cols)) {
      gc <- group_cols

      if (is.null(names(gc))) {
        if (length(gc) < length(fill_levels)) {
          stop(
            "`group_cols` has fewer colors than the number of split groups.",
            call. = FALSE
          )
        }
        gc <- gc[seq_along(fill_levels)]
        names(gc) <- fill_levels
      } else {
        if (!all(fill_levels %in% names(gc))) {
          if (length(gc) >= length(fill_levels)) {
            gc2 <- gc[seq_along(fill_levels)]
            names(gc2) <- fill_levels
            gc <- gc2
          } else {
            stop(
              "Named `group_cols` must contain all split group names.",
              call. = FALSE
            )
          }
        } else {
          gc <- gc[fill_levels]
        }
      }

      return(unname(stats::setNames(gc, fill_levels)))
    }

    fill_scale <- gb$plot$scales$get_scales("fill")
    fill_cols <- NULL

    if (!is.null(fill_scale)) {
      fill_cols <- tryCatch(
        fill_scale$map(fill_levels),
        error = function(e) NULL
      )

      if (is.null(fill_cols) || any(is.na(fill_cols))) {
        fill_cols <- tryCatch(
          fill_scale$palette(length(fill_levels)),
          error = function(e) NULL
        )
      }
    }

    if (is.null(fill_cols) || any(is.na(fill_cols))) {
      fallback_cols <- c("#398fa1", "#f16b00", "#999999", "#CC79A7")
      fill_cols <- rep(fallback_cols, length.out = length(fill_levels))
    }

    names(fill_cols) <- fill_levels
    unname(stats::setNames(fill_cols, fill_levels))
  }

  choose_display_max_break <- function(obs_max_pct, breaks_pct) {
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

  size_from_pct <- function(pct_prop, obs_max_prop,
                            size_range, size_zero, size_power) {
    if (!is.finite(pct_prop) || is.na(pct_prop) || pct_prop <= 0) {
      return(size_zero)
    }

    if (!is.finite(obs_max_prop) || is.na(obs_max_prop) || obs_max_prop <= 0) {
      return(size_zero)
    }

    scaled_prop <- pct_prop / obs_max_prop
    scaled_prop <- max(min(scaled_prop, 1), 0)

    size_range[1] +
      (scaled_prop ^ size_power) *
      (size_range[2] - size_range[1])
  }

  if (!inherits(plot, "ggplot")) {
    stop("`geom_dot_gradient()` must be added to a ggplot object.", call. = FALSE)
  }

  if (!object$direction %in% c(1, -1)) {
    stop("`direction` must be 1 or -1.", call. = FALSE)
  }

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

  if (!is.numeric(object$size_range) || length(object$size_range) != 2 ||
      object$size_range[1] < 0 || object$size_range[2] <= object$size_range[1]) {
    stop(
      "`size_range` must be a numeric vector of length 2 with increasing values.",
      call. = FALSE
    )
  }

  if (!is.numeric(object$size_zero) || length(object$size_zero) != 1 ||
      object$size_zero < 0) {
    stop("`size_zero` must be a non-negative number.", call. = FALSE)
  }

  if (!is.numeric(object$size_power) || length(object$size_power) != 1 ||
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

  violin_layer <- find_violin_layer(plot, object$layer)
  gb <- ggplot2::ggplot_build(plot)
  vdat <- gb$data[[violin_layer]]

  raw_df <- plot$data
  if (!is.data.frame(raw_df)) {
    stop(
      "Cannot extract plotting data from the current VlnPlot object.",
      call. = FALSE
    )
  }

  x_var <- aes_name(plot$mapping$x %||% plot$layers[[violin_layer]]$mapping$x)
  y_var <- aes_name(plot$mapping$y %||% plot$layers[[violin_layer]]$mapping$y)
  fill_var <- aes_name(plot$mapping$fill %||% plot$layers[[violin_layer]]$mapping$fill)

  if (is.null(x_var) || is.null(y_var)) {
    stop("Cannot determine x/y variables from the current plot.", call. = FALSE)
  }

  if (!x_var %in% colnames(raw_df) || !y_var %in% colnames(raw_df)) {
    stop(
      "The x/y variables inferred from the plot are not found in plot$data.",
      call. = FALSE
    )
  }

  df <- raw_df[, unique(c(x_var, y_var, fill_var)), drop = FALSE]
  colnames(df)[colnames(df) == x_var] <- ".x"
  colnames(df)[colnames(df) == y_var] <- ".y"

  if (!is.null(fill_var) && fill_var %in% colnames(df)) {
    colnames(df)[colnames(df) == fill_var] <- ".fill"
  } else {
    df$.fill <- "All"
  }

  df <- df[is.finite(df$.y) & !is.na(df$.x) & !is.na(df$.fill), , drop = FALSE]

  if (!is.factor(df$.x)) {
    df$.x <- factor(df$.x, levels = unique(df$.x))
  }
  if (!is.factor(df$.fill)) {
    df$.fill <- factor(df$.fill, levels = unique(df$.fill))
  }

  x_levels <- levels(df$.x)
  fill_levels <- levels(df$.fill)

  fill_offsets <- subgroup_offsets(length(fill_levels), object$dodge_width)
  names(fill_offsets) <- fill_levels

  fill_cols <- resolve_fill_cols(object$group_cols, fill_levels, gb)
  names(fill_cols) <- fill_levels

  agg_tmp <- stats::aggregate(
    .y ~ .x + .fill,
    data = df,
    FUN = function(z) c(
      pct_expr = mean(z > object$expr_cutoff, na.rm = TRUE),
      avg_expr = mean(z, na.rm = TRUE)
    )
  )

  stat_part <- agg_tmp$.y
  if (is.matrix(stat_part)) {
    stat_df <- as.data.frame(stat_part, stringsAsFactors = FALSE)
  } else {
    stat_df <- as.data.frame(do.call(rbind, stat_part), stringsAsFactors = FALSE)
  }

  sum_df <- cbind(
    agg_tmp[, c(".x", ".fill"), drop = FALSE],
    stat_df
  )

  avg_rng <- range(sum_df$avg_expr, na.rm = TRUE)
  avg_den <- diff(avg_rng)
  if (!is.finite(avg_den) || avg_den <= 0) {
    avg_den <- 1
  }

  q_low_val <- avg_rng[1] + object$low_quantile * avg_den
  q_mid_val <- avg_rng[1] + object$mid_quantile * avg_den
  q_high_val <- avg_rng[1] + object$high_quantile * avg_den

  if (object$direction == 1) {
    alpha_low_end <- object$alpha_min
    alpha_mid <- 0.5
    alpha_high_end <- object$alpha_max
  } else {
    alpha_low_end <- object$alpha_max
    alpha_mid <- 0.5
    alpha_high_end <- object$alpha_min
  }

  sum_df$alpha_val <- vapply(
    sum_df$avg_expr,
    function(y) {
      interp_alpha_piecewise(
        y = y,
        q_low_y = q_low_val,
        q_mid_y = q_mid_val,
        q_high_y = q_high_val,
        alpha_low_end = alpha_low_end,
        alpha_mid = alpha_mid,
        alpha_high_end = alpha_high_end
      )
    },
    numeric(1)
  )

  obs_max_pct <- max(sum_df$pct_expr, na.rm = TRUE) * 100
  if (!is.finite(obs_max_pct) || is.na(obs_max_pct)) {
    obs_max_pct <- 0
  }

  display_max_break <- choose_display_max_break(
    obs_max_pct = obs_max_pct,
    breaks_pct = object$guide_size_breaks
  )

  obs_max_prop <- max(sum_df$pct_expr, na.rm = TRUE)
  if (!is.finite(obs_max_prop) || is.na(obs_max_prop)) {
    obs_max_prop <- 0
  }

  sum_df$size_val <- vapply(
    sum_df$pct_expr,
    function(p) {
      size_from_pct(
        pct_prop = p,
        obs_max_prop = obs_max_prop,
        size_range = object$size_range,
        size_zero = object$size_zero,
        size_power = object$size_power
      )
    },
    numeric(1)
  )

  sum_df$fill_col <- unname(fill_cols[as.character(sum_df$.fill)])
  sum_df$border_col <- sum_df$fill_col

  y_rng <- range(vdat$y, na.rm = TRUE)
  y_span <- diff(y_rng)
  if (!is.finite(y_span) || y_span <= 0) {
    y_span <- 1
  }

  sum_df$y_plot <- y_rng[1] - y_span * object$y_offset
  sum_df$x_num <- as.numeric(sum_df$.x) + unname(fill_offsets[as.character(sum_df$.fill)])

  lower_lim <- min(sum_df$y_plot, na.rm = TRUE) - y_span * 0.12

  anchor_df <- data.frame(
    .x = factor(x_levels[1], levels = x_levels),
    .y_anchor = lower_lim,
    stringsAsFactors = FALSE
  )

  plot_main <- plot +
    ggplot2::geom_blank(
      data = anchor_df,
      mapping = ggplot2::aes(x = .x, y = .y_anchor),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = sum_df,
      mapping = ggplot2::aes(
        x = x_num,
        y = y_plot,
        size = I(size_val),
        fill = I(fill_col),
        colour = I(border_col),
        alpha = I(alpha_val)
      ),
      inherit.aes = FALSE,
      shape = object$dot_shape,
      stroke = object$dot_stroke,
      show.legend = FALSE
    ) +
    ggplot2::coord_cartesian(clip = "off")

  if (!is.null(object$sig_test_all)) {
    df$.group_key <- interaction(df$.x, df$.fill, drop = TRUE, lex.order = TRUE)
    grp_levels <- levels(df$.group_key)

    test_df <- data.frame(
      group_key = grp_levels,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    )

    for (i in seq_along(grp_levels)) {
      g <- grp_levels[i]
      target <- df$.y[df$.group_key == g]
      background <- df$.y[df$.group_key != g]

      p <- tryCatch(
        {
          if (object$sig_test_all == "wilcox") {
            suppressWarnings(
              stats::wilcox.test(
                target,
                background,
                alternative = "greater",
                exact = FALSE
              )$p.value
            )
          } else {
            suppressWarnings(
              stats::t.test(
                target,
                background,
                alternative = "greater"
              )$p.value
            )
          }
        },
        error = function(e) NA_real_
      )

      test_df$p_value[i] <- p
    }

    test_df$p_adj <- stats::p.adjust(test_df$p_value, method = object$sig_adjust)
    test_df$star <- vapply(
      test_df$p_adj,
      p_to_star,
      character(1),
      cutoffs = object$sig_cutoffs
    )

    key_map <- unique(data.frame(
      group_key = interaction(sum_df$.x, sum_df$.fill, drop = TRUE, lex.order = TRUE),
      .x = sum_df$.x,
      .fill = sum_df$.fill,
      x_num = sum_df$x_num,
      y_plot = sum_df$y_plot,
      stringsAsFactors = FALSE
    ))

    star_df <- merge(test_df, key_map, by = "group_key", all.x = TRUE)
    star_df <- star_df[star_df$star != "" & !is.na(star_df$x_num), , drop = FALSE]

    if (nrow(star_df) > 0) {
      star_df$y_star <- star_df$y_plot + y_span * object$star_offset

      plot_main <- plot_main +
        ggplot2::geom_text(
          data = star_df,
          mapping = ggplot2::aes(
            x = x_num,
            y = y_star,
            label = star
          ),
          inherit.aes = FALSE,
          colour = "black",
          size = object$star_size,
          fontface = "bold",
          vjust = 0.5,
          show.legend = FALSE
        )
    }
  }

  y_scale <- plot_main$scales$get_scales("y")
  if (!is.null(y_scale)) {
    y_scale$limits <- c(lower_lim, NA)
    y_scale$oob <- scales::oob_keep
  } else {
    plot_main <- plot_main +
      ggplot2::scale_y_continuous(
        limits = c(lower_lim, NA),
        oob = scales::oob_keep
      )
  }

  if (!isTRUE(object$guide)) {
    return(plot_main)
  }

  guide_breaks_all <- sort(unique(as.numeric(object$guide_size_breaks)))
  guide_breaks_show <- guide_breaks_all[guide_breaks_all <= display_max_break]
  if (!0 %in% guide_breaks_show) {
    guide_breaks_show <- c(0, guide_breaks_show)
  }
  guide_breaks_show <- sort(unique(guide_breaks_show))

  guide_max_break <- max(guide_breaks_show[guide_breaks_show > 0], na.rm = TRUE)
  if (!is.finite(guide_max_break)) {
    guide_max_break <- 100
  }

  size_vals <- vapply(
    guide_breaks_show / 100,
    function(p) {
      if (p <= 0) {
        object$size_zero
      } else {
        scaled_prop <- (p * 100) / guide_max_break
        object$size_range[1] +
          (scaled_prop ^ object$size_power) *
          (object$size_range[2] - object$size_range[1])
      }
    },
    numeric(1)
  )

  size_df <- data.frame(
    x = rep(1, length(guide_breaks_show)),
    y = rev(seq_along(guide_breaks_show)) + 1,
    s = rev(unname(size_vals)),
    lab = paste0(rev(unname(guide_breaks_show)), "%"),
    stringsAsFactors = FALSE
  )

  guide_bar_df <- data.frame(
    x = rep(2, object$guide_bar_n),
    y = seq(1, 7, length.out = object$guide_bar_n),
    val = seq(avg_rng[1], avg_rng[2], length.out = object$guide_bar_n),
    stringsAsFactors = FALSE
  )

  tick_df <- data.frame(
    x = c(2.30, 2.30, 2.30),
    y = c(1, 4, 7),
    lab = sprintf("%.2f", c(avg_rng[1], mean(avg_rng), avg_rng[2])),
    stringsAsFactors = FALSE
  )

  guide_plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 8.3,
      label = "Expressed\npercentage",
      hjust = 0,
      size = object$guide_title_size / 3
    ) +
    ggplot2::geom_point(
      data = size_df,
      mapping = ggplot2::aes(x = x, y = y, size = I(s)),
      shape = 16,
      colour = "black"
    ) +
    ggplot2::geom_text(
      data = size_df,
      mapping = ggplot2::aes(x = x + 0.28, y = y, label = lab),
      hjust = 0,
      size = object$guide_label_size / 3
    ) +
    ggplot2::annotate(
      "text",
      x = 1.85,
      y = 8.3,
      label = "Average\nexpression",
      hjust = 0,
      size = object$guide_title_size / 3
    ) +
    ggplot2::geom_tile(
      data = guide_bar_df,
      mapping = ggplot2::aes(x = x, y = y, fill = val),
      width = 0.20,
      height = 6 / object$guide_bar_n
    ) +
    ggplot2::geom_text(
      data = tick_df,
      mapping = ggplot2::aes(x = x, y = y, label = lab),
      hjust = 0,
      size = object$guide_label_size / 3
    ) +
    ggplot2::scale_fill_gradient(
      low = object$guide_bar_low,
      high = object$guide_bar_high
    ) +
    ggplot2::xlim(0.8, 3.0) +
    ggplot2::ylim(0.5, 8.8) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  patchwork::wrap_plots(
    plot_main,
    guide_plot,
    ncol = 2,
    widths = c(1, object$guide_width)
  )
}