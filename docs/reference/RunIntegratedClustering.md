# Cluster cells on an integrated embedding

Builds a k-nearest-neighbour graph on a corrected embedding and detects
communities for downstream trajectory and pseudobulk workflows.

## Usage

``` r
RunIntegratedClustering(
  sce,
  reduction = NULL,
  cluster_col = "integrated_clusters",
  dims = NULL,
  k = 20,
  algorithm = c("louvain", "walktrap", "leiden"),
  resolution = 0.8,
  seed = 1
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- reduction:

  Integrated reduced-dimension name. If `NULL`, the function prefers a
  stored `integrated_*` embedding.

- cluster_col:

  Column name used to store cluster labels in `colData(sce)`.

- dims:

  Optional dimensions to retain from the embedding.

- k:

  Number of nearest neighbours used to build the graph.

- algorithm:

  Graph clustering algorithm: `"louvain"` or `"walktrap"`.

- resolution:

  Cluster granularity parameter. This is analogous to Seurat's
  `resolution`: higher values typically produce more clusters. It is
  used when the selected graph clustering algorithm supports it.

- seed:

  Random seed used when graph construction requires a fallback.

## Value

A `SingleCellExperiment` with cluster labels added to `colData(sce)`.
