# Run pseudobulk spline or time-series differential expression with edgeR

Fit spline-based edgeR models to identify genes that change over time
within each condition, then optionally return fitted expression curves
for visualization.

## Usage

``` r
run_edger_spline_de(
  pb,
  time_col = "timepoint",
  condition_col = NULL,
  covariates = NULL,
  assay_name = "counts",
  normalization = c("TMM", "upperquartile", "none"),
  df = 3,
  filter_by_expr = TRUE,
  robust = TRUE,
  abundance_trend = FALSE,
  fdr_cutoff = 0.05,
  return_curve_for = c("significant", "all"),
  curve_grid_length = 200
)
```

## Arguments

- pb:

  A pseudobulk `SingleCellExperiment`.

- time_col:

  Sample metadata column containing numeric or ordered time.

- condition_col:

  Optional sample metadata column defining conditions to be analysed
  separately, for example treatment.

- covariates:

  Optional additive covariates to include in each condition-specific
  model.

- assay_name:

  Assay containing pseudobulk counts.

- normalization:

  Library-size normalization method. Choices are:

  - `"TMM"`: edgeR TMM normalization.

  - `"upperquartile"`: edgeR upper-quartile normalization.

  - `"none"`: no normalization factors beyond raw library sizes.

- df:

  Degrees of freedom for the natural spline basis.

- filter_by_expr:

  Logical indicating whether to run
  [`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html).

- robust:

  Logical passed to
  [`edgeR::glmQLFit()`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html).

- abundance_trend:

  Logical passed to
  [`edgeR::glmQLFit()`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html).

- fdr_cutoff:

  FDR cutoff used to define significant genes.

- return_curve_for:

  Which genes should receive fitted spline curves. Choices are:

  - `"significant"`: return curves only for significant genes.

  - `"all"`: return curves for all tested genes.

- curve_grid_length:

  Number of points used when drawing fitted spline curves.

## Value

A list containing condition-specific edgeR fits, tables of temporal DE
genes, and fitted curve data for plotting.
