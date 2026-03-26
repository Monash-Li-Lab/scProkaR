# Plot any stored reduced dimension

Creates a quick scatter plot for an existing reduced dimension such as
an original Seurat UMAP, an integrated latent space, or a UMAP computed
by
[`RunIntegratedUMAP()`](https://monash-li-lab.github.io/scProka/reference/RunIntegratedUMAP.md).

## Usage

``` r
PlotReduction(
  sce,
  reduction,
  colour_by = NULL,
  point_size = 0.4,
  point_alpha = 0.8,
  palette = NULL,
  facet_by = NULL,
  shuffle = FALSE,
  seed = 1
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- reduction:

  Reduced-dimension name to plot.

- colour_by:

  Optional `colData(sce)` column used for colouring points.

- point_size:

  Point size.

- point_alpha:

  Point alpha.

- palette:

  Optional colour palette. For discrete variables, this should be a
  character vector of colours; for continuous variables, it should be a
  gradient vector.

- facet_by:

  Optional `colData(sce)` column used to facet the plot.

- shuffle:

  Whether to shuffle plotting order before drawing points.

- seed:

  Random seed used when `shuffle = TRUE`.

## Value

A `ggplot2` object.
