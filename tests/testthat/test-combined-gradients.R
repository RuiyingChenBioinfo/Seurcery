test_that("dot summaries are calculated separately for each facet", {
  d <- gradient_test_data()
  d$value[d$gene == "small" & d$cluster == "C1"] <- 0
  p <- ggplot2::ggplot(d, ggplot2::aes(cluster, value, fill = gene)) +
    ggplot2::geom_violin() + ggplot2::facet_wrap(~gene, scales = "free_y")
  out <- suppressWarnings(p + geom_vln_gradient(bin = 12) +
                           geom_dot_gradient(guide = FALSE))
  summary <- dot_summary_data(out)
  expect_equal(nrow(summary), 4L)
  expect_equal(sort(summary$avg_expr),
               sort(as.numeric(tapply(d$value, interaction(d$cluster, d$gene), mean))))
  expect_equal(sort(summary$pct_expr), c(0, 100, 100, 100))
  built <- suppressWarnings(ggplot2::ggplot_build(out))
  point_layers <- which(vapply(out$layers, function(l) inherits(l$geom, "GeomPoint"),
                               logical(1)))
  expect_length(point_layers, 1)
  expect_equal(as.integer(table(built$data[[point_layers]]$PANEL)), c(2L, 2L))
})

test_that("both addition orders work when dot guides create a patchwork", {
  p <- gradient_test_plot()
  first_violin <- p + geom_vln_gradient(bin = 12) + geom_dot_gradient(guide = TRUE)
  first_dot <- p + geom_dot_gradient(guide = TRUE) + geom_vln_gradient(bin = 12)
  expect_s3_class(first_violin, "patchwork")
  expect_s3_class(first_dot, "patchwork")
  expect_equal(dot_summary_data(first_violin), dot_summary_data(first_dot))
  original_alpha <- lapply(panel_records(p + geom_vln_gradient(bin = 12)),
                           function(x) vapply(x, `[[`, numeric(1), "alpha"))
  for (out in list(first_violin, first_dot)) {
    leaves <- violin_leaves(out)
    expect_length(leaves, 1)
    expect_true(inherits(leaves[[1]]$layers[[1]]$geom, "GeomSeurceryViolinGradient"))
    combined_alpha <- lapply(panel_records(leaves[[1]]),
                             function(x) vapply(x, `[[`, numeric(1), "alpha"))
    expect_equal(combined_alpha, original_alpha)
    expect_s3_class(render_gradient_plot(out), "gtable")
  }
})

test_that("dots preserve existing size and alpha scales and accept later fill colours", {
  d <- subset(gradient_test_data(), gene == "small")
  p <- ggplot2::ggplot(d, ggplot2::aes(cluster, value, fill = cluster)) +
    ggplot2::geom_violin() +
    ggplot2::geom_point(ggplot2::aes(size = value, alpha = value)) +
    ggplot2::scale_size_continuous(range = c(0.3, 2)) +
    ggplot2::scale_alpha_continuous(range = c(0.2, 0.5))
  old_points <- ggplot2::ggplot_build(p)$data[[2]]
  out <- p + geom_vln_gradient(bin = 10) + geom_dot_gradient(guide = FALSE) +
    ggplot2::scale_fill_manual(values = c(C1 = "#3FA55E", C2 = "#EA8350"))
  built <- ggplot2::ggplot_build(out)
  expect_equal(built$data[[2]]$alpha, old_points$alpha)
  expect_equal(built$data[[2]]$size, old_points$size)
  dot_layer <- which(vapply(out$layers, function(l) {
    is.data.frame(l$data) && all(c("avg_expr", "pct_expr") %in% names(l$data))
  }, logical(1)))
  expect_length(dot_layer, 1)
  expect_setequal(built$data[[dot_layer]]$fill, c("#3FA55E", "#EA8350"))
})

test_that("nested patchworks retain plot order, layout and annotations", {
  one <- function(title) {
    ggplot2::ggplot(subset(gradient_test_data(), gene == "small"),
                    ggplot2::aes(cluster, value, fill = cluster)) +
      ggplot2::geom_violin() + ggplot2::labs(title = title)
  }
  nested <- (one("first") / one("second")) | one("third")
  nested <- nested + patchwork::plot_layout(widths = c(2, 1), guides = "collect") +
    patchwork::plot_annotation(title = "Keep this title", tag_levels = "A")
  changed <- nested + geom_vln_gradient(bin = 10)
  expect_equal(vapply(patchwork_leaves(changed), function(x) x$labels$title,
                      character(1)), c("first", "second", "third"))
  expect_equal(changed$patches$layout, nested$patches$layout)
  expect_equal(changed$patches$annotation, nested$patches$annotation)
  expect_equal(changed[[1]]$patches$layout, nested[[1]]$patches$layout)
  expect_s3_class(render_gradient_plot(changed), "gtable")
})

test_that("dot summaries respect horizontal violin aesthetics and coord_flip", {
  d <- subset(gradient_test_data(), gene == "small")
  horizontal <- ggplot2::ggplot(d, ggplot2::aes(value, cluster, fill = cluster)) +
    ggplot2::geom_violin(orientation = "y")
  flipped <- ggplot2::ggplot(d, ggplot2::aes(cluster, value, fill = cluster)) +
    ggplot2::geom_violin() + ggplot2::coord_flip()
  for (p in list(horizontal, flipped)) {
    out <- p + geom_vln_gradient(bin = 10) + geom_dot_gradient(guide = FALSE)
    summary <- dot_summary_data(out)
    expect_equal(sort(summary$avg_expr), sort(as.numeric(tapply(d$value, d$cluster, mean))))
    expect_s3_class(render_gradient_plot(out), "gtable")
  }
  expect_s3_class((flipped + geom_dot_gradient(guide = FALSE))$coordinates, "CoordFlip")
})
