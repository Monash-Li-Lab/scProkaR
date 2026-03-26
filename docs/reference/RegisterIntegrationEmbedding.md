# Register an external embedding for benchmarking and downstream analysis

Adds a precomputed embedding to `reducedDims(sce)` and records it in
`metadata(sce)$SCProkaR$integration`, allowing custom integration
outputs to be benchmarked alongside native SCProkaR methods.

## Usage

``` r
RegisterIntegrationEmbedding(
  sce,
  embedding,
  method_name,
  reduction_name = NULL,
  metadata = list()
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- embedding:

  A matrix-like cell embedding with one row per cell.

- method_name:

  Name used when recording the embedding.

- reduction_name:

  Optional reducedDim name. Defaults to
  `paste0("integrated_", method_name)`.

- metadata:

  Optional named list describing the external method.

## Value

A `SingleCellExperiment` with the external embedding registered.
