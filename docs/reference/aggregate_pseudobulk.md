# Aggregate single-cell data into pseudobulk samples

Summarize single-cell counts into pseudobulk samples defined by one or
more metadata columns. This is a common first step before sample-level
differential expression analysis with tools such as edgeR.

## Usage

``` r
aggregate_pseudobulk(
  sce,
  sample_cols,
  assay_name = "counts",
  aggregation = c("sum", "mean"),
  min_cells = 1L,
  add_logcounts = TRUE,
  prior_count = 2,
  sample_prefix = "pb_"
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- sample_cols:

  Character vector of `colData(sce)` columns that together define each
  pseudobulk sample, for example
  `c("sample_id", "timepoint", "treatment", "cluster")`.

- assay_name:

  Assay containing the input counts.

- aggregation:

  Aggregation method. Choices are:

  - `"sum"`: sum counts across cells in each pseudobulk sample. This is
    the recommended default for count-based edgeR workflows.

  - `"mean"`: average expression across cells in each pseudobulk sample.
    This can be useful for exploratory summaries but is not ideal for
    count-based NB models.

- min_cells:

  Minimum number of cells required for a pseudobulk sample to be
  retained.

- add_logcounts:

  Logical indicating whether to add a `logcounts` assay using
  `edgeR::cpm(..., log = TRUE)`.

- prior_count:

  Prior count passed to
  [`edgeR::cpm()`](https://rdrr.io/pkg/edgeR/man/cpm.html) when
  `add_logcounts = TRUE`.

- sample_prefix:

  Prefix for the generated pseudobulk column names.

## Value

A pseudobulk `SingleCellExperiment` with one column per aggregated
sample. The output `colData` includes the grouping metadata, `ncells`,
and `lib.size`.
