# Build the TATA cluster graph

Combine topology and temporal consistency into final TATA edge weights,
assign edge direction, and prune weak edges to obtain the cluster-level
TATA graph.

## Usage

``` r
build_tata_graph(
  topology_table,
  time_table,
  clusters,
  embedding,
  timepoint,
  alpha = 1,
  beta = 1,
  direction_threshold = 0.1,
  prune_threshold = 0.01,
  time_weight_mode = c("directional_confidence", "ordered")
)
```

## Arguments

- topology_table:

  Output from
  [`compute_topology_weights()`](https://monash-li-lab.github.io/scProka/reference/compute_topology_weights.md).

- time_table:

  Output from
  [`compute_time_weights()`](https://monash-li-lab.github.io/scProka/reference/compute_time_weights.md).

- clusters:

  Cluster labels for cells.

- embedding:

  Low-dimensional embedding used to compute cluster centroids.

- timepoint:

  Experimental time labels.

- alpha:

  Exponent for topology weights.

- beta:

  Exponent for time weights.

- direction_threshold:

  Threshold for assigning edge direction.

- prune_threshold:

  Minimum TATA weight to retain an edge.

- time_weight_mode:

  How to incorporate time into the final edge weight. Choices are:

  - `"directional_confidence"`: reward strong directional evidence
    regardless of whether the supported direction is `a -> b` or
    `b -> a`.

  - `"ordered"`: use the original ordered-pair formulation
    `(flow_score + 1) / 2`, which is stricter but depends on the
    ordering of `cluster_a` and `cluster_b`.

## Value

A list containing the directed cluster graph, the cluster edge table,
and cluster vertex metadata.
