# Measure the weight of anti-causal edges in a directed graph

Compare each directed graph edge with the median experimental time of
its source and target clusters, and quantify how much total edge weight
points backward in time.

## Usage

``` r
anticausal_edge_mass(
  cluster_graph,
  clusters,
  timepoint,
  tolerance = 0,
  weight_attr = "weight"
)
```

## Arguments

- cluster_graph:

  Directed cluster-level graph as an `igraph`.

- clusters:

  Cluster label for each cell.

- timepoint:

  Experimental time labels for the same cells.

- tolerance:

  A non-negative tolerance for small median-time differences. Directed
  edges with time differences in `[-tolerance, tolerance]` are treated
  as time-neutral rather than anti-causal.

- weight_attr:

  Edge attribute to use as the edge weight.

## Value

A list with:

- `summary`: a one-row `data.frame` containing the total edge weight,
  the anti-causal edge mass, its proportion, and a higher-is-better
  causal edge score.

- `edge_table`: a `data.frame` describing each directed edge.
