# Seurcery 0.8.7

* Add per-panel violin opacity gradients for stacked Seurat plots, including
  `VlnPlot(..., stack = TRUE, flip = TRUE)` with a gene-specific manual fill scale.
* Draw gradient slices after ggplot2 scale training to retain source colors and
  allow fill scales to be added before or after the violin gradient.
* Detect the expression axis for horizontal and vertical violins and keep
  gradient ranges and dot summaries separate across facets.
* Support both addition orders of `geom_vln_gradient()` and
  `geom_dot_gradient()`, including dot-guide panels.
* Traverse nested patchworks while preserving layouts, annotations, spacers,
  and non-violin panels.
* Share layer detection and gradient validation; report invalid settings with
  clearer errors.
* Clarify that the `*_quantile` arguments are range fractions rather than
  empirical quantiles, and correct the documented `min_avg_pct_diff` default
  to 0.2 (20 percentage points).
* Add regression coverage and a reproducible 15-gene Seurat example.
