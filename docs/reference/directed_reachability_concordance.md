# Quantify how well a directed cluster graph follows real time

For ordered cluster pairs with increasing median experimental time,
measure whether the directed cluster graph supports reachability in the
forward direction.

## Usage

``` r
directed_reachability_concordance(
  cluster_graph,
  clusters,
  timepoint,
  min_time_gap = 0,
  weight_by = c("cells", "pairs")
)
```

## Arguments

- cluster_graph:

  Directed cluster-level graph as an `igraph`.

- clusters:

  Cluster label for each cell.

- timepoint:

  Experimental time labels for the same cells.

- min_time_gap:

  Minimum difference in median cluster time required before a pair is
  considered a meaningful forward-time comparison.

- weight_by:

  How to weight cluster-pair comparisons. Choices are:

  - `"cells"`: weight each ordered pair by the product of the two
    cluster sizes.

  - `"pairs"`: weight every ordered pair equally.

## Value

A list with:

- `summary`: a one-row `data.frame` containing the weighted forward
  reachability concordance and the number of ordered cluster pairs.

- `pair_table`: a `data.frame` describing each ordered cluster pair.
