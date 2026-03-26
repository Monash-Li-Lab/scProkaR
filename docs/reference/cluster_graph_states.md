# Cluster graph-connected cell states

Cluster the cell graph into coarse states using a graph-based community
detection method.

## Usage

``` r
cluster_graph_states(
  cell_graph,
  method = c("louvain", "leiden", "walktrap"),
  seed = 1
)
```

## Arguments

- cell_graph:

  An undirected cell-level `igraph`.

- method:

  Clustering method. Choices are:

  - `"louvain"`: fast modularity-based community detection.

  - `"leiden"`: Leiden clustering if supported by the installed
    `igraph`.

  - `"walktrap"`: random-walk clustering that can produce broader
    partitions.

- seed:

  Random seed.

## Value

A factor of inferred cluster labels.

## Examples

``` r
sce <- simulate_tata_sce(n_cells = 400, n_features = 60, seed = 1)
#> Warning: The divergent simulation is designed for larger datasets; consider at least 1000 cells for clearer branch structure.
knn <- build_knn_graph(sce, dimred = "PCA", k = 10)
clusters <- cluster_graph_states(knn$graph, method = "louvain", seed = 1)
table(clusters)
#> clusters
#>  C1  C2  C3  C4  C5  C6 
#>  70  95  65  39 102  29 
```
