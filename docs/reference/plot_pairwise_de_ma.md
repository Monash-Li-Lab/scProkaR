# Plot an MA plot for a pairwise edgeR contrast

Visualize one pairwise contrast from
[`run_edger_pairwise_de()`](https://monash-li-lab.github.io/scProka/reference/run_edger_pairwise_de.md)
as an MA plot using average abundance and log-fold change.

## Usage

``` r
plot_pairwise_de_ma(
  de_result,
  contrast = NULL,
  fdr_cutoff = 0.05,
  point_size = 1.5,
  point_alpha = 0.75
)
```

## Arguments

- de_result:

  Result list returned by
  [`run_edger_pairwise_de()`](https://monash-li-lab.github.io/scProka/reference/run_edger_pairwise_de.md).

- contrast:

  Name of the contrast to plot. If `NULL`, the first available contrast
  is used.

- fdr_cutoff:

  FDR threshold used to define significant genes.

- point_size:

  Point size.

- point_alpha:

  Point alpha.

## Value

A `ggplot2` object.
