test_that("dot summaries warn before omitting missing groups or fills", {
  dat <- data.frame(
    group = rep(c("A", "B"), each = 20),
    fill = rep(c("A", "B"), each = 20),
    value = rep(seq_len(20), 2)
  )
  for (column in c("group", "fill")) {
    missing <- dat
    missing[[column]][1:4] <- NA_character_
    p <- ggplot2::ggplot(missing, ggplot2::aes(group, value, fill = fill)) +
      ggplot2::geom_violin()
    expect_warning(out <- p + geom_dot_gradient(guide = FALSE),
                   "Omitting 4 observations.*dot summaries")
    expect_s3_class(render_gradient_plot(out), "gtable")
  }
})
