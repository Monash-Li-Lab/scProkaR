# Example bacterial single-cell dataset

`sce1` is a packaged `SingleCellExperiment` example dataset for
demonstrating the core `SCProkaR` workflow. It contains one observation
per bacterium and includes three treatment groups (`control`, `PMB0.5`,
and `PMB2`) measured across multiple time points.

## Usage

``` r
sce1
```

## Format

A `SingleCellExperiment` with a counts assay and cell metadata columns:

- sample:

  original sample identifier

- clusters:

  pre-existing cluster label

- treatments:

  treatment group

- timepoints:

  experimental sampling time

## Source

Internal example dataset distributed with `SCProkaR`.
