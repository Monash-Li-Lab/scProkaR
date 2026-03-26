# Simulate a multi-drug branching TATA dataset

Simulate a `SingleCellExperiment` with one shared origin trajectory that
can diverge into three drug-associated fates. The three drug-associated
branches emerge at different real treatment times, making the dataset
useful for checking whether TATA uses time information rather than
feature similarity alone when orienting cluster-level transitions.

## Usage

``` r
simulate_tata_multidrug_sce(
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

The simulated design uses treatment times in hours:

- `drug_A` branch emerges early, around 24 hours

- `drug_C` branch emerges around 48 hours

- `drug_B` branch emerges late, around 96 hours

The geometry is intentionally center-rooted rather than tip-rooted:
`drug_A` and `drug_B` extend in opposite directions from a shared middle
state, and `drug_C` diverges along a third arm. This is meant to
stress-test whether time labels can recover direction when the root is
not visually obvious from topology alone.

## Examples

``` r
sce <- simulate_tata_multidrug_sce(n_cells = 1000, n_features = 80, seed = 1)
table(sce$condition, sce$timepoint)
#>         
#>           0 24 48 72 96 120
#>   drug_A 53 55 54 60 54  55
#>   drug_B 67 49 67 55 67  46
#>   drug_C 62 41 64 46 64  41
```
