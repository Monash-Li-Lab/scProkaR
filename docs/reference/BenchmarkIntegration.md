# Benchmark integration quality with scIB-inspired metrics

Computes a compact R-native metric panel across available reduced
dimensions and summarizes the results with scores and optional plots.

## Usage

``` r
BenchmarkIntegration(
  sce,
  batch_col,
  label_col = NULL,
  methods = NULL,
  metrics = c("pcr_batch", "ilisi", "clisi", "kbet_like", "graph_connectivity", "ari",
    "nmi", "silhouette"),
  return_plots = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- batch_col:

  Column in `colData(sce)` containing batch labels.

- label_col:

  Optional biological label column for conservation metrics.

- methods:

  Optional vector of integration method names or reducedDim names to
  benchmark. When `NULL`, all stored integration results are used.

- metrics:

  Metric names to calculate.

- return_plots:

  If `TRUE`, return summary ggplot objects.

  Any embedding already present in `reducedDims(sce)` can be benchmarked
  by passing its reduced-dimension name in `methods`, and external
  embeddings can be registered with
  [`RegisterIntegrationEmbedding()`](https://monash-li-lab.github.io/scProka/reference/RegisterIntegrationEmbedding.md).

## Value

A list with `scores`, `ranking`, and optional `plots`.
