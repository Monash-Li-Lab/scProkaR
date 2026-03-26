# Simulate a divergent TATA dataset

Simulate a `SingleCellExperiment` containing single-bacterium
observations on a divergent two-way trajectory. Time `0` is centred at
the origin of the latent trajectory, and later timepoints move away from
that central root into two branches. This makes the simulated dataset
useful for demonstrating TATA on a divergent trajectory where the
earliest state appears in the middle of a UMAP.

## Usage

``` r
simulate_tata_sce(n_cells = 5000, n_features = 100, n_pcs = 20, seed = 1)

simulate_tata_divergent_sce(
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

The returned object contains:

- `counts` and `logcounts`

- `colData` with `cell_id`, `timepoint`, `condition`, `simulated_state`,
  and `simulated_branch`

- `reducedDims(sce)$PCA`

- `reducedDims(sce)$UMAP`

- `reducedDims(sce)$truth`

## Examples

``` r
sce <- simulate_tata_sce(n_cells = 500, n_features = 60, seed = 1)
#> Warning: The divergent simulation is designed for larger datasets; consider at least 1000 cells for clearer branch structure.
dim(sce)
#> [1]  60 500
```
