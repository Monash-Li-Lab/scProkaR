# Run a full SCProkaR differential expression workflow

Aggregate single-cell counts into pseudobulk samples and then dispatch
to a pairwise or time-series edgeR workflow.

## Usage

``` r
run_scproka_de(
  sce,
  mode = c("pairwise", "time_series"),
  sample_cols,
  assay_name = "counts",
  aggregation = c("sum", "mean"),
  min_cells = 10L,
  normalization = c("TMM", "upperquartile", "none"),
  ...
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- mode:

  Workflow mode. Choices are:

  - `"pairwise"`: pseudobulk followed by pairwise edgeR contrasts.

  - `"time_series"`: pseudobulk followed by spline-based time-series
    edgeR.

- sample_cols:

  Metadata columns used to define pseudobulk samples.

- assay_name:

  Assay containing input counts.

- aggregation:

  Pseudobulk aggregation method. Choices are:

  - `"sum"`: sum counts within each pseudobulk sample.

  - `"mean"`: average counts within each pseudobulk sample.

- min_cells:

  Minimum number of cells required per pseudobulk sample.

- normalization:

  Pseudobulk normalization method. Choices are:

  - `"TMM"`: edgeR TMM normalization.

  - `"upperquartile"`: edgeR upper-quartile normalization.

  - `"none"`: no normalization factors beyond raw library sizes.

- ...:

  Additional arguments passed to
  [`run_edger_pairwise_de()`](https://monash-li-lab.github.io/scProka/reference/run_edger_pairwise_de.md)
  or
  [`run_edger_spline_de()`](https://monash-li-lab.github.io/scProka/reference/run_edger_spline_de.md),
  depending on `mode`.

## Value

A list with the pseudobulk object, the DE result, and workflow settings.
