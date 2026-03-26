# Compute local temporal order metrics on a cell graph

Summarize how often pseudotime agrees with the known direction of time
on local cell-cell graph edges.

## Usage

``` r
local_temporal_order_metrics(pseudotime, adjacency, timepoint)
```

## Arguments

- pseudotime:

  Numeric vector of inferred pseudotime values.

- adjacency:

  Sparse adjacency matrix for the cell-cell graph.

- timepoint:

  Experimental time labels for the same cells.

## Value

A one-row `data.frame` with the number of compared edges, local time
concordance, local temporal inversion rate, a higher-is-better temporal
direction score, and the fraction of edges tied in pseudotime.
