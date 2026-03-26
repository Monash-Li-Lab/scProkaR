# Simulate an extreme spiral-plus-divergent TATA dataset

Simulate a stress-test `SingleCellExperiment` whose `UMAP` embedding
follows a spiral-like manifold that transitions into two divergent
tails. This is useful for visually stress-testing TATA in an
intentionally extreme geometry.

## Usage

``` r
simulate_tata_spiral_sce(
  n_cells = 5000,
  n_features = 100,
  n_pcs = 20,
  seed = 1
)
```

## Arguments

- n_cells:

  Number of cells to simulate.

- n_features:

  Number of features to simulate.

- n_pcs:

  Number of principal components to store.

- seed:

  Random seed.

## Value

A `SingleCellExperiment`.

## Details

The embedding is constructed directly in UMAP-like space and paired with
smooth expression programs so that neighboring cells along the path
remain transcriptionally similar.

## Examples

``` r
sce <- simulate_tata_spiral_sce(n_cells = 1000, n_features = 80, seed = 1)
SingleCellExperiment::reducedDimNames(sce)
#> [1] "PCA"   "UMAP"  "truth"
```
