# Compute time consistency for cluster-cluster edges

For each inter-cluster cell-cell graph edge, compute `sign(t_b - t_a)`
after orienting the cluster pair as `(cluster_a, cluster_b)`. The mean
sign is the temporal flow score.

## Usage

``` r
compute_time_weights(adjacency, clusters, timepoint)
```

## Arguments

- adjacency:

  Sparse cell-cell adjacency matrix.

- clusters:

  Cluster labels for cells.

- timepoint:

  Experimental time labels for cells.

## Value

A data frame of temporal flow scores and time weights.
