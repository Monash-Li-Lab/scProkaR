# Filter pseudobulk samples by size

Filter a pseudobulk `SingleCellExperiment` using minimum numbers of
cells and library size.

## Usage

``` r
filter_pseudobulk_samples(pb, min_cells = 10L, min_lib_size = NULL)
```

## Arguments

- pb:

  A pseudobulk `SingleCellExperiment`, typically returned by
  [`aggregate_pseudobulk()`](https://monash-li-lab.github.io/scProka/reference/aggregate_pseudobulk.md).

- min_cells:

  Minimum number of cells required.

- min_lib_size:

  Minimum library size required.

## Value

A filtered pseudobulk `SingleCellExperiment`.
