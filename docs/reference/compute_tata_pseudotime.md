# Compute TATA pseudotime

Compute a simple cluster-level pseudotime on the directed TATA graph
using shortest-path distances from the root cluster with inverse edge
weights, then propagate that pseudotime to cells using within-cluster
time ranks. If some clusters are not reachable in the directed graph,
TATA falls back to the undirected abstract topology for those clusters
so late branches are not dropped simply because an intermediate
direction is ambiguous.

## Usage

``` r
compute_tata_pseudotime(
  cluster_graph,
  sce,
  cluster_col = "cluster",
  time_col = "timepoint",
  root_cluster = NULL,
  within_cluster_fraction = 0.5,
  eps = 1e-08
)
```

## Arguments

- cluster_graph:

  Directed cluster-level TATA graph.

- sce:

  A `SingleCellExperiment`.

- cluster_col:

  Column in `colData(sce)` containing inferred clusters.

- time_col:

  Column in `colData(sce)` containing experimental time labels.

- root_cluster:

  Optional root cluster. If `NULL`, the root is chosen as the cluster
  with the earliest median timepoint.

- within_cluster_fraction:

  Scale factor used to spread cells within each cluster after assigning
  cluster-level pseudotime.

- eps:

  Small constant to avoid division by zero.

## Value

A list with root cluster, cluster pseudotime table, raw cell pseudotime,
scaled cell pseudotime, and the within-cluster step size.
