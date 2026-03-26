# Plot fitted spline expression curves from time-series DE

Plot fitted spline curves for selected genes, optionally overlaid across
conditions. When a requested gene was significant in one condition but
not another, the non-significant condition is still drawn using a dotted
line.

## Usage

``` r
plot_time_series_deg_curves(
  spline_result,
  genes = NULL,
  top_n = 12L,
  scales = "free_y",
  ncol = NULL
)
```

## Arguments

- spline_result:

  Result list returned by
  [`run_edger_spline_de()`](https://monash-li-lab.github.io/scProka/reference/run_edger_spline_de.md).

- genes:

  Optional character vector of genes to plot. If `NULL`, the top genes
  ranked by FDR are used.

- top_n:

  Number of genes to plot when `genes = NULL`.

- scales:

  Facet scale setting passed to
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).

- ncol:

  Optional number of columns used by
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).

## Value

A `ggplot2` object.
