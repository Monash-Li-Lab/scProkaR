# Plot a volcano plot for a pairwise edgeR contrast

Visualize one pairwise contrast from
[`run_edger_pairwise_de()`](https://monash-li-lab.github.io/scProka/reference/run_edger_pairwise_de.md)
as a volcano plot using log-fold change and FDR.

## Usage

``` r
plot_pairwise_de_volcano(
  de_result,
  contrast = NULL,
  fdr_cutoff = 0.05,
  lfc_cutoff = 0,
  top_n_labels = 10L,
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

- lfc_cutoff:

  Absolute log-fold-change threshold used to highlight larger effects.

- top_n_labels:

  Number of top genes to label by FDR. Set to `0` to omit labels.

- point_size:

  Point size.

- point_alpha:

  Point alpha.

## Value

A `ggplot2` object.
