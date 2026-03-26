# Build a cell-cell kNN graph

Build an undirected cell-cell k-nearest-neighbour graph from a
low-dimensional representation.

## Usage

``` r
build_knn_graph(x, dimred = "PCA", dims = NULL, k = 15)
```

## Arguments

- x:

  A `SingleCellExperiment` or a numeric matrix with cells in rows.

- dimred:

  Name of the reduced dimension to use when `x` is a
  `SingleCellExperiment`.

- dims:

  Optional subset of dimensions to use.

- k:

  Number of neighbours.

## Value

A list containing an `igraph` object, a sparse adjacency matrix, the
embedding used, neighbour indices, neighbour distances, and kNN
metadata.

## Examples

``` r
sce <- simulate_tata_sce(n_cells = 400, n_features = 60, seed = 1)
#> Warning: The divergent simulation is designed for larger datasets; consider at least 1000 cells for clearer branch structure.
knn <- build_knn_graph(sce, dimred = "PCA", k = 10)
igraph::vcount(knn$graph)
#> [1] 400
```
