# Seurcery

<img src="Public/Logo-Seurcery.png" align="right" width="115" height="135">

An R package for extending Seurat visualizations. Version 0.8.7 adds per-gene
gradients to stacked violin plots and improves combined violin and dot gradients.

The [v0.8.6 tutorial](Public/Seurcery_v0.8.6_tutorial.html) describes the
previous release. See the examples below and [NEWS.md](NEWS.md) for v0.8.7.

<details> <summary> Note about package development </summary>

Seurcery is actively being developed. You may occasionally encounter bugs or minor documentation issues. GitHub issues are welcome for bug reports, usage questions, and feature suggestions.

</details>

## Installation

Install the v0.8.7 source archive after installing its dependencies (`grid` is
included with R):

```r
install.packages(c("ggplot2", "patchwork", "rlang"))
install.packages("Seurcery-Seurcery_v0.8.7.tar.gz", repos = NULL, type = "source")
```

## Stacked violin gradients

The following uses 15 genes from your Seurat object. Name the palette with `gl`
so that each color remains associated with its gene if the plotting order changes.

```r
library(ggplot2)
library(Seurcery)

# rds: your Seurat object; gl: 15 unique gene names present in rds
stopifnot(length(gl) == 15L, !anyDuplicated(gl))
gene_cols <- setNames(
  c(rep("#3FA55E", 2), rep("#50B5D8", 5), rep("#EA8350", 8)),
  gl
)

p <- Seurat::VlnPlot(
  rds,
  features = gl,
  group.by = "seurat_clusters",
  pt.size = 0,
  stack = TRUE,
  flip = TRUE
) +
  scale_fill_manual(values = gene_cols)

p + geom_vln_gradient(alpha_min = 0.05, alpha_max = 1)
```

Each gene panel has its own expression range for the opacity gradient. Low
expression is lighter and high expression is darker, while the fill color comes
from the original plot's scale. Use `direction = -1` to reverse the gradient.
Both `flip = TRUE` and `flip = FALSE` are supported, as are horizontal ggplot
violins and `coord_flip()`.

Violin gradients are drawn after scale training, so `scale_fill_manual()` may
also be added after `geom_vln_gradient()`. For a patchwork of separate feature
plots, use patchwork's `&` operator to apply a fill scale to every child plot.
Existing patchwork layouts, annotations, spacers, and non-violin panels are kept.

## Combining violin and dot gradients

Add optional dots showing the expressed percentage (size) and relative average
expression (opacity) for each group. Summaries and opacity ranges are calculated
separately for each gene panel. For exact expressed percentages at
`expr_cutoff = 0`, set `add.noise = FALSE` in `Seurat::VlnPlot()` when creating
`p`; Seurat's default plotting noise can make zero values slightly positive.

```r
p +
  geom_vln_gradient(alpha_min = 0.05, alpha_max = 1) +
  geom_dot_gradient(guide = FALSE)

# The reverse order is also supported.
p +
  geom_dot_gradient(guide = FALSE) +
  geom_vln_gradient(alpha_min = 0.05, alpha_max = 1)
```

Set `guide = TRUE` (the default) to include the dot guide. Both addition orders
also work with this guide: adding a violin gradient visits the violin panels
and leaves guide-only panels intact. Dots follow the source fill scale unless
`group_cols` is supplied. For horizontal violins, dot summaries are placed to
the left of the expression range.

The existing argument names `low_quantile`, `mid_quantile`, and `high_quantile`
denote fractions of the expression range, not empirical quantiles. Violin
anchors use the panel's plotted density range; dot anchors use the range of
group average expression within that panel. For example, `mid_quantile = 0.5`
means the midpoint of the relevant range.

For a reproducible example using Seurat's `pbmc_small`, run:

```sh
Rscript examples/stacked-gradients.R /path/to/stacked-gradients.pdf
```

The script plots 15 genes with both flip settings and both gradient orders,
with and without dot guides.

## Citation
Chen, R. (2026). Seurcery: An R package designed for further analysis and visualization based on utilities of R package Seurat (Version v0.8.7) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22953630

## Contact
* Ruiying Chen (chenruiying21@mails.ucas.ac.cn)
