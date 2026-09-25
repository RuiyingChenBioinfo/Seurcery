gradient_test_data <- function() {
  d <- expand.grid(
    cell = seq_len(40), cluster = c("C1", "C2"),
    gene = c("small", "large"), KEEP.OUT.ATTRS = FALSE
  )
  d$value <- (d$cell / 40 + as.integer(d$cluster) / 10) *
    ifelse(d$gene == "small", 1, 100)
  d
}

gradient_test_plot <- function() {
  ggplot2::ggplot(gradient_test_data(),
                  ggplot2::aes(cluster, value, fill = gene)) +
    ggplot2::geom_violin() +
    ggplot2::facet_wrap(~gene, scales = "free_y")
}

walk_grobs <- function(g) {
  descendants <- c(if (!is.null(g$children)) as.list(g$children),
                   if (!is.null(g$grobs)) g$grobs)
  c(list(g), unlist(lapply(descendants, walk_grobs), recursive = FALSE))
}

# Record the actual rendered polygon vertices and RGBA values. Checking grobs
# catches panel recycling and frozen colours that ggplot_build alone misses.
polygon_records <- function(g) {
  polys <- Filter(function(x) inherits(x, "polygon"), walk_grobs(g))
  unname(unlist(lapply(polys, function(p) {
    fills <- p$gp$fill
    if (is.null(fills) || all(is.na(fills))) return(list())
    ids <- p$id
    if (is.null(ids) && !is.null(p$id.lengths)) {
      ids <- rep(seq_along(p$id.lengths), p$id.lengths)
    }
    if (is.null(ids)) ids <- rep(1L, length(p$x))
    groups <- unique(ids)
    lapply(seq_along(groups), function(i) {
      fill <- fills[(i - 1L) %% length(fills) + 1L]
      if (is.na(fill)) return(NULL)
      rgba <- grDevices::col2rgb(fill, alpha = TRUE)[, 1]
      alpha <- rgba[4] / 255
      if (!is.null(p$gp$alpha)) alpha <- alpha * p$gp$alpha[1]
      list(x = as.numeric(p$x[ids == groups[i]]),
           y = as.numeric(p$y[ids == groups[i]]),
           rgb = toupper(grDevices::rgb(rgba[1], rgba[2], rgba[3],
                                        maxColorValue = 255)),
           alpha = unname(alpha))
    })
  }), recursive = FALSE))
}

gradient_pbmc_small <- function() {
  env <- new.env(parent = emptyenv())
  utils::data("pbmc_small", package = "SeuratObject", envir = env)
  env$pbmc_small
}

panel_records <- function(plot) {
  gt <- render_gradient_plot(plot)
  indices <- which(grepl("^panel($|-)", gt$layout$name))
  lapply(gt$grobs[indices], polygon_records)
}

patchwork_leaves <- function(p) {
  if (!inherits(p, "patchwork")) return(list(p))
  unlist(lapply(seq_len(length(p)), function(i) patchwork_leaves(p[[i]])),
         recursive = FALSE)
}

violin_leaves <- function(p) {
  Filter(function(x) any(vapply(x$layers, function(l) {
    any(grepl("Violin", class(l$geom)))
  }, logical(1))), patchwork_leaves(p))
}

dot_summary_data <- function(p) {
  pieces <- unlist(lapply(violin_leaves(p), function(leaf) {
    lapply(leaf$layers, function(l) {
      if (is.data.frame(l$data) &&
          all(c("avg_expr", "pct_expr") %in% names(l$data))) l$data else NULL
    })
  }), recursive = FALSE)
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(data.frame())
  do.call(rbind, pieces)
}

render_gradient_plot <- function(p) {
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output, width = 7, height = 7)
  device <- grDevices::dev.cur()
  on.exit({
    grDevices::dev.off(device)
    unlink(output)
  }, add = TRUE)
  if (inherits(p, "patchwork")) patchwork::patchworkGrob(p)
  else ggplot2::ggplotGrob(p)
}
