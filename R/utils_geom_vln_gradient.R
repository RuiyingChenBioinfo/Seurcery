#' Add gradient fill polygons to violin plots
#'
#' Creates a lightweight layer-like object that can be added to a ggplot object
#' with `+`. The object is handled by the S3 method
#' `ggplot_add.geom_vln_gradient()`, which reconstructs violin polygons from an
#' existing violin layer and fills them with an alpha gradient along the y axis.
#'
#' @param bin Number of vertical bins used to reconstruct the gradient fill.
#'   Must be a single number greater than or equal to 2.
#' @param alpha_max Maximum alpha value.
#' @param alpha_min Minimum alpha value.
#' @param direction Gradient direction. Use `1` for low to high alpha from low
#'   y to high y, and `-1` for the reverse.
#' @param layer Optional index of the violin layer in the ggplot object. If
#'   `NULL`, the function tries to detect the first violin-like layer
#'   automatically.
#' @param outline Logical. Whether to redraw violin outlines.
#' @param outline_size Optional outline width. If `NULL`, uses the linewidth in
#'   the source violin layer when available.
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
geom_vln_gradient <- function(bin = 50,
                              alpha_max = 1,
                              alpha_min = 0,
                              direction = 1,
                              layer = NULL,
                              outline = TRUE,
                              outline_size = NULL,
                              min_width_frac = 0.02,
                              high_quantile = 1,
                              mid_quantile = 0.5,
                              low_quantile = 0) {
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
#' @param plot A ggplot object.
#' @param object_name Not used.
#'
#' @return A ggplot object with reconstructed violin polygons.
#' @export
ggplot_add.geom_vln_gradient <- function(object, plot, object_name) {
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      y
    } else {
      x
    }
  }

  if (!inherits(plot, "ggplot")) {
    stop("`geom_vln_gradient()` must be added to a ggplot object.", call. = FALSE)
  }

  if (!is.numeric(object$bin) || length(object$bin) != 1 || object$bin < 2) {
    stop("`bin` must be a single number >= 2.", call. = FALSE)
  }

  if (!is.numeric(object$alpha_min) || !is.numeric(object$alpha_max) ||
      object$alpha_min < 0 || object$alpha_max > 1 ||
      object$alpha_min > object$alpha_max) {
    stop("Need 0 <= alpha_min <= alpha_max <= 1.", call. = FALSE)
  }

  if (!object$direction %in% c(1, -1)) {
    stop("`direction` must be 1 or -1.", call. = FALSE)
  }

  if (!is.numeric(object$min_width_frac) ||
      length(object$min_width_frac) != 1 ||
      object$min_width_frac < 0 ||
      object$min_width_frac >= 1) {
    stop("`min_width_frac` must be in [0, 1).", call. = FALSE)
  }

  q_vals <- c(object$low_quantile, object$mid_quantile, object$high_quantile)
  if (any(!is.numeric(q_vals)) || any(q_vals < 0) || any(q_vals > 1)) {
    stop(
      "`low_quantile`, `mid_quantile`, and `high_quantile` must all be in [0, 1].",
      call. = FALSE
    )
  }

  if (!(object$low_quantile <= object$mid_quantile &&
        object$mid_quantile <= object$high_quantile)) {
    stop("Need low_quantile <= mid_quantile <= high_quantile.", call. = FALSE)
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

  prep_curve <- function(y, x) {
    df <- data.frame(y = y, x = x)
    df <- df[is.finite(df$y) & is.finite(df$x), , drop = FALSE]
    if (nrow(df) == 0) {
      return(NULL)
    }

    df <- stats::aggregate(x ~ y, data = df, FUN = mean)
    df <- df[order(df$y), , drop = FALSE]
    if (nrow(df) < 2) {
      return(NULL)
    }

    df
  }

  build_full_outline <- function(left_df, right_df, col, id) {
    data.frame(
      x = c(left_df$x, rev(right_df$x), left_df$x[1]),
      y = c(left_df$y, rev(right_df$y), left_df$y[1]),
      outline_id = id,
      col = col,
      stringsAsFactors = FALSE
    )
  }

  build_split_outline <- function(side_df, x0, col, id) {
    data.frame(
      x = c(rep(x0, nrow(side_df)), rev(side_df$x), x0),
      y = c(side_df$y, rev(side_df$y), side_df$y[1]),
      outline_id = id,
      col = col,
      stringsAsFactors = FALSE
    )
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

  violin_layer <- find_violin_layer(plot, object$layer)
  gb <- ggplot2::ggplot_build(plot)
  vdat <- gb$data[[violin_layer]]

  need_cols <- c("x", "y", "xmin", "xmax", "violinwidth", "group")
  if (!all(need_cols %in% colnames(vdat))) {
    stop(
      "The selected layer does not contain the columns needed to reconstruct violin shapes.",
      call. = FALSE
    )
  }

  if (!"PANEL" %in% colnames(vdat)) {
    vdat$PANEL <- 1
  }
  if (!"fill" %in% colnames(vdat)) {
    vdat$fill <- "grey70"
  }
  if (!"colour" %in% colnames(vdat) && !"color" %in% colnames(vdat)) {
    vdat$colour <- "black"
  }
  if (!"colour" %in% colnames(vdat) && "color" %in% colnames(vdat)) {
    vdat$colour <- vdat$color
  }

  y_rng <- range(vdat$y, na.rm = TRUE)
  y_den <- y_rng[2] - y_rng[1]
  if (!is.finite(y_den) || y_den <= 0) {
    y_den <- 1
  }

  breaks <- seq(y_rng[1], y_rng[2], length.out = object$bin + 1)
  mids <- (breaks[-1] + breaks[-length(breaks)]) / 2

  q_low_y <- y_rng[1] + object$low_quantile * y_den
  q_mid_y <- y_rng[1] + object$mid_quantile * y_den
  q_high_y <- y_rng[1] + object$high_quantile * y_den

  if (object$direction == 1) {
    alpha_low_end <- object$alpha_min
    alpha_mid <- 0.5
    alpha_high_end <- object$alpha_max
  } else {
    alpha_low_end <- object$alpha_max
    alpha_mid <- 0.5
    alpha_high_end <- object$alpha_min
  }

  alpha_vals <- vapply(
    mids,
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

  alpha_vals[alpha_vals < 0] <- 0
  alpha_vals[alpha_vals > 1] <- 1

  vdat$.xkey <- paste(vdat$PANEL, sprintf("%.10f", vdat$x), sep = "__")
  x_group_n <- tapply(vdat$group, vdat$.xkey, function(z) length(unique(z)))

  split_one <- split(vdat, interaction(vdat$PANEL, vdat$group, drop = TRUE))

  poly_list <- list()
  outline_list <- list()
  pid <- 1
  oid <- 1

  for (gdat in split_one) {
    gdat <- gdat[is.finite(gdat$y) & is.finite(gdat$violinwidth), , drop = FALSE]
    if (nrow(gdat) < 2) {
      next
    }

    gdat <- gdat[order(gdat$y), , drop = FALSE]

    x0 <- unique(gdat$x)[1]
    xkey <- unique(gdat$.xkey)[1]
    n_same_x <- x_group_n[[xkey]]

    x_left <- gdat$x - gdat$violinwidth * (gdat$x - gdat$xmin)
    x_right <- gdat$x + gdat$violinwidth * (gdat$xmax - gdat$x)

    left_df <- prep_curve(gdat$y, x_left)
    right_df <- prep_curve(gdat$y, x_right)
    if (is.null(left_df) || is.null(right_df)) {
      next
    }

    fill_col <- unique(stats::na.omit(gdat$fill))[1] %||% "grey70"
    line_col <- unique(stats::na.omit(gdat$colour))[1] %||% "black"

    lw <- object$outline_size %||%
      if ("linewidth" %in% colnames(gdat)) unique(gdat$linewidth)[1] else NULL %||%
      0.3

    if (n_same_x == 1) {
      yr <- range(gdat$y, na.rm = TRUE)
      full_width_max <- max(right_df$x - left_df$x, na.rm = TRUE)
      width_cut <- full_width_max * object$min_width_frac

      for (i in seq_len(object$bin)) {
        y0 <- breaks[i]
        y1 <- breaks[i + 1]
        ys0 <- max(y0, yr[1])
        ys1 <- min(y1, yr[2])

        if (ys1 <= ys0) {
          next
        }

        xl0 <- stats::approx(left_df$y, left_df$x, xout = ys0, rule = 2, ties = mean)$y
        xl1 <- stats::approx(left_df$y, left_df$x, xout = ys1, rule = 2, ties = mean)$y
        xr0 <- stats::approx(right_df$y, right_df$x, xout = ys0, rule = 2, ties = mean)$y
        xr1 <- stats::approx(right_df$y, right_df$x, xout = ys1, rule = 2, ties = mean)$y

        local_width <- max(c(xr0 - xl0, xr1 - xl1), na.rm = TRUE)
        if (!is.finite(local_width) || local_width <= width_cut) {
          next
        }

        poly_list[[pid]] <- data.frame(
          x = c(xl0, xl1, xr1, xr0),
          y = c(ys0, ys1, ys1, ys0),
          poly_id = paste0("poly_", pid),
          fill_col = fill_col,
          alpha_val = alpha_vals[i],
          stringsAsFactors = FALSE
        )
        pid <- pid + 1
      }

      if (isTRUE(object$outline)) {
        outline_list[[oid]] <- build_full_outline(
          left_df = left_df,
          right_df = right_df,
          col = line_col,
          id = paste0("outline_", oid)
        )
        outline_list[[oid]]$lw <- lw
        oid <- oid + 1
      }

    } else if (n_same_x == 2) {
      grp_id <- unique(gdat$group)[1]
      side_df <- if (grp_id %% 2 == 1) left_df else right_df
      yr <- range(gdat$y, na.rm = TRUE)

      half_width_max <- max(abs(side_df$x - x0), na.rm = TRUE)
      width_cut <- half_width_max * object$min_width_frac

      for (i in seq_len(object$bin)) {
        y0 <- breaks[i]
        y1 <- breaks[i + 1]
        ys0 <- max(y0, yr[1])
        ys1 <- min(y1, yr[2])

        if (ys1 <= ys0) {
          next
        }

        xs0 <- stats::approx(side_df$y, side_df$x, xout = ys0, rule = 2, ties = mean)$y
        xs1 <- stats::approx(side_df$y, side_df$x, xout = ys1, rule = 2, ties = mean)$y

        local_width <- max(c(abs(xs0 - x0), abs(xs1 - x0)), na.rm = TRUE)
        if (!is.finite(local_width) || local_width <= width_cut) {
          next
        }

        poly_list[[pid]] <- data.frame(
          x = c(x0, x0, xs1, xs0),
          y = c(ys0, ys1, ys1, ys0),
          poly_id = paste0("poly_", pid),
          fill_col = fill_col,
          alpha_val = alpha_vals[i],
          stringsAsFactors = FALSE
        )
        pid <- pid + 1
      }

      if (isTRUE(object$outline)) {
        outline_list[[oid]] <- build_split_outline(
          side_df = side_df,
          x0 = x0,
          col = line_col,
          id = paste0("outline_", oid)
        )
        outline_list[[oid]]$lw <- lw
        oid <- oid + 1
      }

    } else {
      stop(
        "Current version supports regular violin or two-group split violin only.",
        call. = FALSE
      )
    }
  }

  if (length(poly_list) == 0) {
    stop("Failed to reconstruct violin polygons.", call. = FALSE)
  }

  poly_df <- do.call(rbind, poly_list)

  plot$layers[[violin_layer]]$aes_params$fill <- NA
  plot$layers[[violin_layer]]$aes_params$colour <- NA
  plot$layers[[violin_layer]]$aes_params$color <- NA

  plot <- plot +
    ggplot2::geom_polygon(
      data = poly_df,
      mapping = ggplot2::aes(
        x = x,
        y = y,
        group = poly_id,
        fill = I(fill_col),
        alpha = I(alpha_val)
      ),
      inherit.aes = FALSE,
      colour = NA,
      show.legend = FALSE
    )

  if (isTRUE(object$outline) && length(outline_list) > 0) {
    outline_df <- do.call(rbind, outline_list)

    plot <- plot +
      ggplot2::geom_path(
        data = outline_df,
        mapping = ggplot2::aes(
          x = x,
          y = y,
          group = outline_id,
          colour = I(col),
          linewidth = I(lw)
        ),
        inherit.aes = FALSE,
        lineend = "round",
        show.legend = FALSE
      )
  }

  plot
}