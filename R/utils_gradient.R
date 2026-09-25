# Shared validation and traversal for both gradient additions.
.seurcery_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

.seurcery_validate_alpha <- function(object) {
  if (!.seurcery_scalar(object$alpha_min) ||
      !.seurcery_scalar(object$alpha_max) || object$alpha_min < 0 ||
      object$alpha_max > 1 || object$alpha_min > object$alpha_max) {
    stop("Need 0 <= alpha_min <= alpha_max <= 1 (finite scalars).", call. = FALSE)
  }
  anchors <- list(object$low_quantile, object$mid_quantile, object$high_quantile)
  if (!all(vapply(anchors, .seurcery_scalar, logical(1)))) {
    stop("Gradient anchors must be finite numeric scalars.", call. = FALSE)
  }
  anchors <- unlist(anchors)
  if (any(anchors < 0 | anchors > 1) || any(diff(anchors) < 0)) {
    stop("Need 0 <= low_quantile <= mid_quantile <= high_quantile <= 1.", call. = FALSE)
  }
  invisible(TRUE)
}

.seurcery_find_violin_layer <- function(plot, layer = NULL, required = TRUE) {
  hits <- which(vapply(plot$layers, function(x) {
    any(grepl("Violin", class(x$geom), ignore.case = TRUE))
  }, logical(1)))
  if (!is.null(layer)) {
    if (!.seurcery_scalar(layer) || layer != floor(layer) ||
        !layer %in% hits) {
      stop("`layer` must be the index of an existing violin layer.", call. = FALSE)
    }
    return(as.integer(layer))
  }
  if (length(hits)) return(hits[1L])
  if (required) stop("Cannot find a violin layer in this plot.", call. = FALSE)
  NULL
}

.seurcery_apply_gradient <- function(plot, object, add_fun) {
  n_changed <- 0L
  walk <- function(p) {
    if (inherits(p, "patchwork")) {
      # Public indexing preserves order, nested layouts, annotations and spacers.
      for (i in seq_len(length(p))) p[[i]] <- walk(p[[i]])
      return(p)
    }
    if (!inherits(p, "ggplot")) return(p)
    if (is.null(.seurcery_find_violin_layer(p, required = FALSE))) return(p)
    n_changed <<- n_changed + 1L
    add_fun(p, object)
  }
  out <- walk(plot)
  if (!n_changed) {
    stop("The plot must contain at least one violin layer.", call. = FALSE)
  }
  out
}

# Work in canonical coordinates (x = group, y = expression), then swap back.
.seurcery_flip_data <- function(data, flipped) {
  if (!isTRUE(flipped)) return(data)
  old <- c("x", "y", "xmin", "ymin", "xmax", "ymax", "xend", "yend")
  new <- c("y", "x", "ymin", "xmin", "ymax", "xmax", "yend", "xend")
  idx <- match(names(data), old)
  names(data)[!is.na(idx)] <- new[idx[!is.na(idx)]]
  data
}
