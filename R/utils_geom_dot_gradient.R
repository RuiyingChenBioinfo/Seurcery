#' Add dot summaries with alpha gradients below violin plots
#'
#' Creates a lightweight layer-like object that can be added to a ggplot object
#' with `+`. The corresponding S3 method `ggplot_add.geom_dot_gradient()`
#' extracts plotting data from an existing violin plot, calculates expressed
#' percentage and average expression for each x and split group, and draws
#' summary dots with gradient alpha, optional significance labels, and an
#' optional side guide panel.
#'
#' This function supports both a single ggplot violin plot and a patchwork
#' object returned by `Seurat::VlnPlot()` when multiple features are combined.
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
geom_dot_gradient <- function(
  expr_cutoff = 0,
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
  guide_bar_high = "grey15"
) {
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

#' Add dot summaries with alpha gradients below violin plots
#'
#' S3 method for adding a `"geom_dot_gradient"` object to a ggplot.
#'
#' @param object An object created by [geom_dot_gradient()].
#' @param plot A ggplot object or a patchwork object.
#' @param object_name Not used.
#'
#' @return A ggplot object or a patchwork object with summary dots.
#' @importFrom ggplot2 ggplot_add
#' @method ggplot_add geom_dot_gradient
#' @export
ggplot_add.geom_dot_gradient <- function(object, plot, object_name) {
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
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

  is_patchwork_like <- function(p) {
    inherits(p, "patchwork") && !is.null(p$patches)
  }

  apply_to_patchwork <- function(pw, add_fun) {
    strip_patchwork <- function(x) {
      y <- x
      class(y) <- setdiff(class(y), "patchwork")
      y$patches <- NULL
      y
    }

    main_plot <- strip_patchwork(pw)
    all_plots <- c(list(main_plot), pw$patches$plots)
    all_plots <- lapply(all_plots, add_fun)

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

  add_single_plot <- function(plot_single) {
    find_violin_layer <- function(plot_obj, layer = NULL) {
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

    interp_alpha_piecewise <- function(
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

    size_from_pct <- function(
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

    if (!inherits(plot_single, "ggplot")) {
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
      stop("`size_zero` must be a non-negative number.", call. = FALSE)
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

    violin_layer <- find_violin_layer(plot_single, object$layer)
    gb <- ggplot2::ggplot_build(plot_single)
    vdat <- gb$data[[violin_layer]]

    raw_df <- plot_single$data
    if (!is.data.frame(raw_df)) {
      stop(
        "Cannot extract plotting data from the current VlnPlot object.",
        call. = FALSE
      )
    }

    x_var <- aes_name(plot_single$mapping$x %||% plot_single$layers[[violin_layer]]$mapping$x)
    y_var <- aes_name(plot_single$mapping$y %||% plot_single$layers[[violin_layer]]$mapping$y)
    fill_var <- aes_name(plot_single$mapping$fill %||% plot_single$layers[[violin_layer]]$mapping$fill)

    if (is.null(x_var) || is.null(y_var)) {
      stop("Cannot determine x/y variables from the current plot.", call. = FALSE)
    }

    if (!x_var %in% colnames(raw_df) || !y_var %in% colnames(raw_df)) {
      stop(
        "The x/y variables inferred from the plot are not found in plot$data.",
        call. = FALSE
      )
    }

    if (!is.null(fill_var) && !fill_var %in% colnames(raw_df)) {
      fill_var <- NULL
    }

    dat <- raw_df[, unique(c(x_var, y_var, fill_var)), drop = FALSE]
    colnames(dat)[colnames(dat) == x_var] <- ".x"
    colnames(dat)[colnames(dat) == y_var] <- ".y"

    if (!is.null(fill_var)) {
      colnames(dat)[colnames(dat) == fill_var] <- ".fill"
      dat$.fill <- as.character(dat$.fill)
      fill_levels <- if (is.factor(raw_df[[fill_var]])) {
        levels(raw_df[[fill_var]])
      } else {
        unique(as.character(dat$.fill))
      }
      fill_levels <- fill_levels[fill_levels %in% unique(dat$.fill)]
    } else {
      dat$.fill <- "All"
      fill_levels <- "All"
    }

    x_levels <- if (is.factor(raw_df[[x_var]])) {
      levels(raw_df[[x_var]])
    } else {
      unique(as.character(raw_df[[x_var]]))
    }
    dat$.x <- as.character(dat$.x)
    x_levels <- x_levels[x_levels %in% unique(dat$.x)]

    dat <- dat[is.finite(dat$.y) & !is.na(dat$.x) & !is.na(dat$.fill), , drop = FALSE]
    if (nrow(dat) == 0) {
      stop("No valid data found for dot summary.", call. = FALSE)
    }

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
    display_max_pct <- choose_display_max_break(obs_max_pct, object$guide_size_breaks)
    display_max_prop <- display_max_pct / 100

    sum_df$size_val <- vapply(
      sum_df$pct_expr / 100,
      function(p) {
        size_from_pct(
          pct_prop = p,
          obs_max_prop = display_max_prop,
          size_range = object$size_range,
          size_zero = object$size_zero,
          size_power = object$size_power
        )
      },
      numeric(1)
    )

    fill_cols <- resolve_fill_cols(object$group_cols, fill_levels, gb)
    names(fill_cols) <- fill_levels

    y_rng <- range(vdat$y, na.rm = TRUE)
    y_span <- y_rng[2] - y_rng[1]
    if (!is.finite(y_span) || y_span <= 0) {
      y_span <- 1
    }

    x_num <- seq_along(x_levels)
    names(x_num) <- x_levels

    offs <- subgroup_offsets(length(fill_levels), object$dodge_width)
    names(offs) <- fill_levels

    sum_df$x_plot <- x_num[as.character(sum_df$.x)] + offs[as.character(sum_df$.fill)]
    sum_df$y_plot <- y_rng[1] - object$y_offset * y_span
    sum_df$fill_col <- fill_cols[as.character(sum_df$.fill)]

    star_df <- NULL

    if (!is.null(object$sig_test_all) && length(fill_levels) > 1) {
      p_tab <- list()

      for (xx in x_levels) {
        sub_x <- dat[dat$.x == xx, , drop = FALSE]
        if (nrow(sub_x) == 0) {
          next
        }

        for (ff in fill_levels) {
          g1 <- sub_x$.y[sub_x$.fill == ff]
          g2 <- sub_x$.y[sub_x$.fill != ff]

          if (length(g1) == 0 || length(g2) == 0) {
            next
          }

          p_val <- tryCatch(
            {
              if (object$sig_test_all == "wilcox") {
                stats::wilcox.test(g1, g2)$p.value
              } else {
                stats::t.test(g1, g2)$p.value
              }
            },
            error = function(e) NA_real_
          )

          p_tab[[length(p_tab) + 1]] <- data.frame(
            .x = xx,
            .fill = ff,
            p = p_val,
            stringsAsFactors = FALSE
          )
        }
      }

      if (length(p_tab) > 0) {
        p_df <- do.call(rbind, p_tab)
        p_df$p_adj <- stats::p.adjust(p_df$p, method = object$sig_adjust)
        p_df$label <- vapply(p_df$p_adj, p_to_star, character(1), cutoffs = object$sig_cutoffs)
        p_df <- p_df[nzchar(p_df$label), , drop = FALSE]

        if (nrow(p_df) > 0) {
          p_df$x_plot <- x_num[p_df$.x] + offs[p_df$.fill]
          p_df$y_plot <- y_rng[1] - (object$y_offset - object$star_offset) * y_span
          star_df <- p_df
        }
      }
    }

    bottom_limit <- y_rng[1] - (object$y_offset + 0.08) * y_span

    main_plot <- plot_single +
      ggplot2::geom_point(
        data = sum_df,
        mapping = ggplot2::aes(
          x = x_plot,
          y = y_plot,
          size = I(size_val),
          fill = I(fill_col),
          alpha = I(alpha_val)
        ),
        inherit.aes = FALSE,
        shape = object$dot_shape,
        stroke = object$dot_stroke,
        colour = "black",
        show.legend = FALSE
      ) +
      ggplot2::expand_limits(y = bottom_limit) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(5.5, 10, 20, 5.5)
      )

    if (!is.null(star_df) && nrow(star_df) > 0) {
      main_plot <- main_plot +
        ggplot2::geom_text(
          data = star_df,
          mapping = ggplot2::aes(
            x = x_plot,
            y = y_plot,
            label = label
          ),
          inherit.aes = FALSE,
          colour = "black",
          size = object$star_size,
          vjust = 0.5,
          fontface = "plain",
          show.legend = FALSE
        )
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
          size_from_pct(
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

    guide_plot_size <- ggplot2::ggplot(guide_size_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(
        ggplot2::aes(size = I(size_val)),
        shape = object$dot_shape,
        stroke = object$dot_stroke,
        fill = "grey70",
        colour = "black"
      ) +
      ggplot2::geom_text(
        data = guide_size_df,
        ggplot2::aes(
          x = 1.55,
          y = y,
          label = paste0(pct, "%")
        ),
        inherit.aes = FALSE,
        hjust = 0,
        size = object$guide_label_size / 3
      ) +
      ggplot2::annotate(
        "text",
        x = 1,
        y = max(guide_size_df$y) + 1.1,
        label = "Expressed percentage",
        size = object$guide_title_size / 3,
        fontface = "bold"
      ) +
      ggplot2::xlim(0.6, 3.1) +
      ggplot2::ylim(0.5, max(guide_size_df$y) + 1.8) +
      ggplot2::theme_void()

    guide_plot_bar <- ggplot2::ggplot(guide_bar_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_tile(
        ggplot2::aes(fill = fill_val),
        height = 1 / object$guide_bar_n,
        width = 0.28
      ) +
      ggplot2::scale_fill_gradientn(
        colours = c(object$guide_bar_low, object$guide_bar_high),
        limits = c(0, 1)
      ) +
      ggplot2::annotate(
        "text",
        x = 1,
        y = 1.12,
        label = "Average expression",
        size = object$guide_title_size / 3,
        fontface = "bold"
      ) +
      ggplot2::annotate(
        "text",
        x = 1.38,
        y = 0,
        label = "Low",
        hjust = 0,
        size = object$guide_label_size / 3
      ) +
      ggplot2::annotate(
        "text",
        x = 1.38,
        y = 1,
        label = "High",
        hjust = 0,
        size = object$guide_label_size / 3
      ) +
      ggplot2::xlim(0.7, 2.1) +
      ggplot2::ylim(-0.05, 1.18) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none")

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

  if (is_patchwork_like(plot)) {
    return(apply_to_patchwork(plot, add_single_plot))
  }

  add_single_plot(plot)
}