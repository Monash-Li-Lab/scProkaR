# Run pairwise pseudobulk differential expression with edgeR

Fit a quasi-likelihood negative-binomial model with edgeR and test
pairwise contrasts between groups such as timepoints, treatments, or
composite treatment-timepoint groups.

## Usage

``` r
run_edger_pairwise_de(
  pb,
  assay_name = "counts",
  group_col = "timepoint",
  covariates = NULL,
  combine_group_cols = NULL,
  contrast_type = c("reference", "all", "manual"),
  reference_level = NULL,
  manual_pairs = NULL,
  normalization = c("TMM", "upperquartile", "none"),
  filter_by_expr = TRUE,
  robust = TRUE,
  abundance_trend = FALSE,
  fdr_cutoff = 0.05
)
```

## Arguments

- pb:

  A pseudobulk `SingleCellExperiment`.

- assay_name:

  Assay containing pseudobulk counts.

- group_col:

  Sample metadata column defining the primary grouping, for example
  `"timepoint"`.

- covariates:

  Optional metadata columns to include as additive covariates in the
  design matrix, for example treatment, batch, or donor.

- combine_group_cols:

  Optional metadata columns to combine into a single composite group via
  [`interaction()`](https://rdrr.io/r/base/interaction.html), for
  example `c("treatment", "timepoint")`.

- contrast_type:

  Contrast construction strategy. Choices are:

  - `"reference"`: compare every non-reference group against one
    reference level. This is the default.

  - `"all"`: test all pairwise group contrasts.

  - `"manual"`: use `manual_pairs` exactly as supplied.

- reference_level:

  Optional reference level used when `contrast_type = "reference"`. By
  default the earliest or first group is used.

- manual_pairs:

  Optional two-column matrix or data frame of manual contrasts with
  `group1` and `group2` values. Each test is `group1 - group2`.

- normalization:

  Library-size normalization method. Choices are:

  - `"TMM"`: edgeR TMM normalization.

  - `"upperquartile"`: edgeR upper-quartile normalization.

  - `"none"`: no normalization factors beyond raw library sizes.

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

  FDR cutoff used to mark significant genes.

## Value

A list containing the fitted edgeR objects, contrast matrix, result
tables, and analysis settings.
