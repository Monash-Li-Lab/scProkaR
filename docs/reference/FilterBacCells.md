# Filter cells using bacterial QC thresholds

Applies simple threshold-based filtering to QC metrics produced by
[`RunBacQC()`](https://monash-li-lab.github.io/scProka/reference/RunBacQC.md).

## Usage

``` r
FilterBacCells(
  sce,
  min_counts = NULL,
  min_features = NULL,
  max_rrna_fraction = NULL,
  custom_filter = NULL
)
```

## Arguments

- sce:

  A `SingleCellExperiment` with QC metrics already present.

- min_counts:

  Minimum allowed total counts.

- min_features:

  Minimum allowed number of detected features.

- max_rrna_fraction:

  Maximum allowed rRNA fraction.

- custom_filter:

  Optional custom filter. May be a logical vector, a column name in
  `colData(sce)`, or a function returning a logical vector.

## Value

A filtered `SingleCellExperiment`.
