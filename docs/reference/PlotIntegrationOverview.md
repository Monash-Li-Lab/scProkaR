# Plot Seurat-style integration result panels

Creates a standard set of integration diagnostic plots inspired by
common Seurat workflows: unintegrated and integrated views coloured by
batch and, optionally, a biological label, plus a split-panel view of
the integrated embedding.

## Usage

``` r
PlotIntegrationOverview(
  sce,
  unintegrated_reduction = "umap",
  integrated_reduction = NULL,
  batch_col,
  label_col = NULL,
  split_by = NULL,
  point_size = 0.4,
  point_alpha = 0.8,
  shuffle = FALSE,
  seed = 1
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- unintegrated_reduction:

  Reduction name for the pre-integration view.

- integrated_reduction:

  Reduction name for the integrated view. If `NULL`, the function
  prefers a stored `umap_integrated_*` reduction and otherwise falls
  back to `integrated_*`.

- batch_col:

  Batch column in `colData(sce)`.

- label_col:

  Optional biological label column in `colData(sce)`.

- split_by:

  Optional metadata column used to facet the integrated plot. This is
  the `SCProkaR` equivalent of Seurat's `split.by`.

- point_size:

  Point size passed to
  [`PlotReduction()`](https://monash-li-lab.github.io/scProka/reference/PlotReduction.md).

- point_alpha:

  Point alpha passed to
  [`PlotReduction()`](https://monash-li-lab.github.io/scProka/reference/PlotReduction.md).

- shuffle:

  Whether to shuffle plotting order.

- seed:

  Random seed used when `shuffle = TRUE`.

## Value

A named list of `ggplot2` objects.
