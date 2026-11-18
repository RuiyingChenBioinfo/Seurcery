# Seurcery

<img src="Public/Logo-Seurcery.png" align="right" width="115" height="135">

An R package for extending Seurat visualizations. Version 0.90 adds
`FeatureLayerPlot()` for displaying multiple genes on aligned thin panels,
alongside the existing violin and dot gradients.

See the [v0.90 tutorial](Public/Seurcery_v0.90_tutorial.html), the examples below,
and [NEWS.md](NEWS.md). The [v0.8.6 tutorial](Public/Seurcery_v0.8.6_tutorial.html)
is retained for the earlier release.

<details> <summary> Note about package development </summary>

Seurcery is actively being developed. You may occasionally encounter bugs or minor documentation issues. GitHub issues are welcome for bug reports, usage questions, and feature suggestions.

</details>

## Installation

Install the v0.90 source archive after installing its dependencies (`grid` is
included with R). `FeatureLayerPlot()` requires SeuratObject >= 5.0.0 and
ggplot2 >= 3.5.0:

```r
install.packages(c("ggplot2", "SeuratObject", "patchwork", "rlang"))
install.packages("Seurcery_0.90.tar.gz", repos = NULL, type = "source")
```

Install Seurat to run the examples below: `install.packages("Seurat")`.

## Feature layers

`FeatureLayerPlot()` shows a metadata reference map on the left and gene
expression on aligned thin panels on the right. Expression is read directly
from the selected assay's normalized `data` layer. All gene panels share one
color scale labeled **Expression level**.

The following example uses the bundled `pbmc_small` data and its stored PCA
coordinates, so no data download or new embedding is needed:

```r
library(Seurat)
library(ggplot2)
library(Seurcery)

data("pbmc_small", package = "SeuratObject")
rds <- pbmc_small
rds$cell_type <- Idents(rds)
gl <- c("LYZ", "MS4A1", "CD3E", "NKG7")

FeatureLayerPlot(
  rds,
  features = gl,
  group.by = "cell_type",
  reduction = "pca",
  assay = "RNA"
)
```

For an object with UMAP coordinates, use `reduction = "umap"` (the default).
The selected reduction must already exist in the object; its first two
dimensions are used. The lower-left marker has two equal, perpendicular black
lines and is labeled automatically, for example **UMAP** or **PCA**.

![FeatureLayerPlot on PBMC3K UMAP coordinates](Public/FeatureLayerPlot_UMAP.png)

Categorical metadata uses Seurat's default discrete colors. Pass a named manual
scale to `group.scale` to update both the reference points and their legend:

```r
pal <- setNames(c("#3FA55E", "#50B5D8", "#EA8350"), levels(rds$cell_type))

FeatureLayerPlot(
  rds, gl, group.by = "cell_type", reduction = "pca",
  group.scale = scale_color_manual(values = pal)
)

# scale_fill_manual(values = pal) is also supported.
FeatureLayerPlot(
  rds, gl, group.by = "nCount_RNA", reduction = "pca"
)
```

Numeric metadata uses an embedded Reds palette and a separate continuous
colorbar. Convert numeric category codes to a factor when they should be
treated as categories. The left title preserves the metadata column name.
`feature.cols` controls the gene gradient, and `text.color` controls text
colors (default `"black"`).

| Argument | Default | Usage |
| --- | --- | --- |
| `order` | `TRUE` | Draw higher values last; `FALSE` keeps the original cell order. Applies to genes and numeric metadata. |
| `ncol` | `2` | Number of gene-panel columns. |
| `gap.x` | `0.18` | Horizontal gap; must be greater than zero. |
| `gap.y` | `0.16` | Vertical gap; must be at least `max(0.10, title.gap + 0.037)`. |
| `title.gap` | `0.04` | Gene-title offset above each panel; must be at least `0.037`. |

Gaps use normalized plot coordinates. Smaller values make the layout more
compact; invalid gaps produce an error before plotting. For example:

```r
FeatureLayerPlot(
  rds, gl, group.by = "cell_type", reduction = "pca",
  order = TRUE, gap.x = 0.06, gap.y = 0.11, title.gap = 0.04
)
```

Missing assay, reduction, normalized data, or requested features produce an
error. Non-finite embedding coordinates also produce an error. Partial missing
or non-finite expression and metadata values generate counted warnings and
remain missing-colored points (grey by default); they are never replaced with
zero. An entirely missing metadata column or requested gene produces an error.
Zero expression is valid. Embedding cells absent from the selected assay are
omitted from all panels with a warning; no shared cells produces an error.

For a reproducible PDF containing categorical, continuous, manual-color, and
compact examples, run after installing Seurat:

```sh
Rscript examples/feature-layers.R /path/to/feature-layers.pdf
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
Chen, R. (2026). Seurcery: An R package designed for further analysis and visualization based on utilities of R package Seurat (Version v0.90) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22953630

## Contact
* Ruiying Chen (chenruiying21@mails.ucas.ac.cn)
