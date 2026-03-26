# Integrate microbial single-cell batches

Runs one of several integration backends and stores the corrected
embedding in `reducedDims(sce)`.

## Usage

``` r
IntegrateBacData(
  sce,
  batch_col,
  method = c("mnn", "harmony"),
  feature_set = c("hvg", "all"),
  dims = 1:30,
  assay_name = "counts",
  integrated_name = NULL,
  ...
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- batch_col:

  Column in `colData(sce)` containing batch assignments.

- method:

  Integration backend: `"mnn"` or `"harmony"`.

- feature_set:

  Either `"hvg"` or `"all"`.

- dims:

  Integer vector of dimensions to retain.

- assay_name:

  Assay used for preprocessing.

- integrated_name:

  Optional name for the output reduced dimension.

- ...:

  Additional backend-specific arguments.

## Value

A `SingleCellExperiment` with an integrated reduced dimension.
