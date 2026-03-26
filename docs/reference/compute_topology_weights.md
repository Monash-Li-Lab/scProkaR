# Compute topology-based cluster connectivity weights

For each cluster pair `(a, b)`, count observed inter-cluster edges and
compare them to a simple degree-aware expectation.

## Usage

``` r
compute_topology_weights(adjacency, clusters, eps = 1e-08)
```

## Arguments

- adjacency:

  Sparse cell-cell adjacency matrix.

- clusters:

  Cluster labels for cells.

- eps:

  Small constant to avoid division by zero.

## Value

A data frame of observed edges, expected edges, and topology weights.
