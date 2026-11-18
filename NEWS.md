# Seurcery 0.9.0

* Add `FeatureLayerPlot()` for normalized gene expression on aligned thin
  panels, with a categorical or numeric metadata reference map.
* Support Seurat discrete colors, manual color or fill scales, an embedded
  Reds palette for numeric metadata, and a shared expression colorbar.
* Mark each labeled break on the expression and continuous metadata colorbars
  in `FeatureLayerPlot()`.
* Select existing embeddings with `reduction` and label them beside equal,
  perpendicular black coordinate markers.
* Add expression ordering, horizontal and vertical panel gaps, and gene-title
  spacing while retaining the established default layout.
* Warn about partial missing expression or metadata values without replacing
  them with zero; report unusable data and invalid spacing with explicit errors.
* Warn with a count when `geom_dot_gradient()` omits observations with
  non-finite expression or missing group or fill values from dot summaries.
* Show zero with the lowest expression colour when every selected value is
  zero, and report ambiguous gene and cell matches across normalized layers.
* Keep dot summaries visible when group means are equal and align dots and
  significance labels with their violins when split groups are missing.
* Remove white seams from the average expression colorbar in `geom_dot_gradient()`.
* Display `seurat_clusters` with its original cluster IDs on the map and legend.
* Use the processed `pbmc3k.final` dataset in the tutorial and show each of
  the three plotting functions in the README.

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
