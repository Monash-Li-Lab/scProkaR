# Run the full TATA workflow

Run Time Aware Trajectory Analysis (TATA) end-to-end on a
`SingleCellExperiment`.

## Usage

``` r
run_tata(
  sce,
  dimred = "PCA",
  dims = NULL,
  time_col = "timepoint",
  cluster_col = "cluster",
  k = 15,
  cluster_method = c("louvain", "leiden", "walktrap"),
  alpha = 1,
  beta = 1,
  direction_threshold = 0.1,
  prune_threshold = 0.01,
  time_weight_mode = c("directional_confidence", "ordered"),
  do_pseudotime = TRUE,
  root_cluster = NULL,
  n_pcs_if_missing = 20,
  seed = 1
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- dimred:

  Reduced dimension to use for graph construction.

- dims:

  Optional subset of dimensions.

- time_col:

  Column containing experimental time labels.

- cluster_col:

  Name of the column where inferred cluster labels will be stored.

- k:

  Number of neighbours for the cell-cell kNN graph.

- cluster_method:

  Graph clustering method. Choices are:

  - `"louvain"`: fast modularity-based clustering and the default
    option.

  - `"leiden"`: Leiden community detection if available in the installed
    `igraph`; otherwise TATA falls back to Louvain.

  - `"walktrap"`: random-walk clustering that can give slightly coarser
    states.

- alpha:

  Exponent for topology weights.

- beta:

  Exponent for time weights.

- direction_threshold:

  Threshold for direction assignment.

- prune_threshold:

  Minimum TATA weight to retain an abstract graph edge.

- time_weight_mode:

  How the time signal contributes to the final edge score. Choices are:

  - `"directional_confidence"`: use direction strength regardless of
    whether the label ordering is `a -> b` or `b -> a`; this is the
    recommended option for multi-branch graphs.

  - `"ordered"`: use the original ordered-pair score
    `(flow_score + 1) / 2`, which is stricter but more sensitive to the
    arbitrary ordering of cluster labels.

- do_pseudotime:

  Logical indicating whether to compute pseudotime.

- root_cluster:

  Optional root cluster for pseudotime.

- n_pcs_if_missing:

  Number of PCs to compute if `dimred` is missing.

- seed:

  Random seed.

## Value

A list containing the updated `sce`, cell graph, adjacency matrix,
cluster graph, cluster edge table, cluster pseudotime, cell pseudotime,
and parameters used.

## Examples

``` r
sce <- simulate_tata_sce(n_cells = 500, n_features = 60, seed = 1)
#> Warning: The divergent simulation is designed for larger datasets; consider at least 1000 cells for clearer branch structure.
tata <- run_tata(sce, dimred = "PCA", k = 12, seed = 1)
names(tata)
#>  [1] "sce"                    "cell_graph"             "adjacency_matrix"      
#>  [4] "cluster_graph"          "cluster_edge_table"     "cluster_vertex_table"  
#>  [7] "cluster_pseudotime"     "cell_pseudotime"        "cell_pseudotime_scaled"
#> [10] "parameters"             "topology_table"         "time_table"            
```
