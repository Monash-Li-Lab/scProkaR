# Normalize pseudobulk samples for exploratory analysis

Compute edgeR normalization factors and add a `logcounts` assay for
plotting and exploratory analysis.

## Usage

``` r
normalize_pseudobulk(
  pb,
  assay_name = "counts",
  method = c("TMM", "upperquartile", "none"),
  output_assay = "logcounts",
  prior_count = 2
)
```

## Arguments

- pb:

  A pseudobulk `SingleCellExperiment`.

- assay_name:

  Assay containing pseudobulk counts.

- method:

  Normalization method. Choices are:

  - `"TMM"`: trimmed mean of M values normalization via
    `edgeR::calcNormFactors(..., method = "TMM")`.

  - `"upperquartile"`: upper-quartile normalization via edgeR.

  - `"none"`: do not apply normalization factors; still compute
    `logcounts` from the raw library sizes.

- output_assay:

  Name of the output log-scale assay.

- prior_count:

  Prior count passed to `edgeR::cpm(..., log = TRUE)`.

## Value

The input pseudobulk object with updated `colData` normalization fields
and a log-scale assay.
