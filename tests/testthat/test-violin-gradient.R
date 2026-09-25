test_that("each gene facet gets a full alpha range and its own manual colour", {
  cols <- c(small = "#3FA55E", large = "#EA8350")
  p <- gradient_test_plot() + geom_vln_gradient(
    bin = 24, alpha_min = 0.6, alpha_max = 0.9, outline = FALSE,
    min_width_frac = 0
  ) + ggplot2::scale_fill_manual(values = cols)
  panels <- panel_records(p)
  expect_length(panels, 2)
  for (i in seq_along(panels)) {
    a <- vapply(panels[[i]], `[[`, numeric(1), "alpha")
    rgb <- unique(vapply(panels[[i]], `[[`, character(1), "rgb"))
    expect_gt(length(a), 20)
    expect_gte(min(a), 0.6 - 1 / 255)
    expect_lte(max(a), 0.9 + 1 / 255)
    expect_lt(min(a), 0.65)
    expect_gt(max(a), 0.85)
    expression_position <- vapply(panels[[i]], function(x) mean(x$y), numeric(1))
    expect_gt(stats::cor(expression_position, a, method = "spearman"), 0.98)
    expect_identical(rgb, unname(cols[i]))
  }
})

test_that("gradient addition does not mutate its source or discard layers", {
  p <- gradient_test_plot() + ggplot2::geom_point(size = 0.3)
  before <- ggplot2::ggplot_build(p)$data
  original_geom <- p$layers[[1]]$geom
  original_params <- p$layers[[1]]$aes_params
  out <- p + geom_vln_gradient(bin = 10)
  expect_identical(p$layers[[1]]$geom, original_geom)
  expect_identical(p$layers[[1]]$aes_params, original_params)
  expect_equal(ggplot2::ggplot_build(p)$data, before)
  expect_length(out$layers, length(p$layers))
  expect_equal(ggplot2::ggplot_build(out)$data[[2]], before[[2]])
  expect_s3_class(render_gradient_plot(out), "gtable")
})

test_that("direction reverses alpha without changing shapes or colours", {
  p <- gradient_test_plot()
  forward <- panel_records(p + geom_vln_gradient(
    bin = 12, outline = FALSE, direction = 1, min_width_frac = 0))
  reverse <- panel_records(p + geom_vln_gradient(
    bin = 12, outline = FALSE, direction = -1, min_width_frac = 0))
  for (i in seq_along(forward)) {
    af <- vapply(forward[[i]], `[[`, numeric(1), "alpha")
    ar <- vapply(reverse[[i]], `[[`, numeric(1), "alpha")
    expect_equal(af + ar, rep(1, length(af)), tolerance = 1 / 255)
    expect_equal(lapply(forward[[i]], `[[`, "x"),
                 lapply(reverse[[i]], `[[`, "x"))
    expect_equal(lapply(forward[[i]], `[[`, "y"),
                 lapply(reverse[[i]], `[[`, "y"))
  }
})

test_that("overlapping ordinary violins retain both halves", {
  d <- expand.grid(value = seq(0.1, 2, length.out = 40),
                   condition = c("A", "B"))
  p <- ggplot2::ggplot(d, ggplot2::aes("C1", value, fill = condition)) +
    ggplot2::geom_violin(position = "identity") +
    ggplot2::scale_fill_manual(values = c(A = "#3FA55E", B = "#EA8350"))
  records <- panel_records(p + geom_vln_gradient(
    bin = 12, outline = FALSE, min_width_frac = 0))[[1]]
  for (col in c("#3FA55E", "#EA8350")) {
    by_colour <- Filter(function(x) identical(x$rgb, col), records)
    expect_gt(length(by_colour), 0)
    # With one discrete x level, the category centre is panel coordinate 0.5.
    for (slice in by_colour) {
      expect_lt(min(slice$x), 0.5)
      expect_gt(max(slice$x), 0.5)
    }
  }
})

test_that("layer data and aesthetics override global data", {
  actual <- gradient_test_data()
  global <- data.frame(unrelated = c("X", "Y"), bad = c(1000, 2000))
  p <- ggplot2::ggplot(global, ggplot2::aes(unrelated, bad)) +
    ggplot2::geom_violin(data = actual,
                        ggplot2::aes(cluster, value, fill = gene),
                        inherit.aes = FALSE) +
    ggplot2::facet_wrap(~gene, scales = "free_y")
  out <- p + geom_vln_gradient(bin = 12)
  expect_equal(ggplot2::ggplot_build(out)$data[[1]],
               ggplot2::ggplot_build(p)$data[[1]])
  expect_length(panel_records(out), 2)
  dots <- out + geom_dot_gradient(guide = FALSE)
  summary <- dot_summary_data(dots)
  expect_equal(nrow(summary), 4L)
  expect_equal(sort(summary$avg_expr),
               sort(as.numeric(tapply(actual$value,
                 interaction(actual$cluster, actual$gene), mean))))
})

test_that("coordinate flips and logarithmic scales are applied only once", {
  p <- gradient_test_plot() + ggplot2::scale_y_log10() + ggplot2::coord_flip()
  out <- p + geom_vln_gradient(bin = 12, outline = FALSE)
  expect_s3_class(out$coordinates, "CoordFlip")
  expect_equal(ggplot2::ggplot_build(out)$data[[1]],
               ggplot2::ggplot_build(p)$data[[1]])
  records <- unlist(panel_records(out), recursive = FALSE)
  expect_gt(length(records), 20)
  expect_true(all(vapply(records, function(r) {
    all(is.finite(c(r$x, r$y))) &&
      all(r$x >= -1e-6 & r$x <= 1 + 1e-6) &&
      all(r$y >= -1e-6 & r$y <= 1 + 1e-6)
  }, logical(1))))
})

test_that("invalid parameters produce clear errors", {
  p <- gradient_test_plot()
  expect_error(p + geom_vln_gradient(bin = NA_real_), "bin")
  expect_error(p + geom_vln_gradient(bin = Inf), "bin")
  expect_error(p + geom_vln_gradient(alpha_min = 0.9, alpha_max = 0.6), "alpha")
  expect_error(p + geom_vln_gradient(direction = 0), "direction")
  expect_error(p + geom_vln_gradient(layer = 99), "layer")
})
