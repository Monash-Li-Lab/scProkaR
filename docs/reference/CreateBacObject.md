# Create a standardized bacterial single-cell object

`CreateBacObject()` standardizes raw microbial single-cell counts into a
`SingleCellExperiment`, aligns cell and feature metadata, and records
package provenance used by downstream SCProkaR workflows.

## Usage

``` r
CreateBacObject(
  x,
  counts_assay = "counts",
  cell_metadata = NULL,
  feature_metadata = NULL,
  sample_col = NULL,
  batch_col = NULL,
  condition_col = NULL,
  time_col = NULL,
  organism = "bacteria",
  seurat_assay = NULL,
  seurat_layer = "counts",
  transfer_reductions = TRUE
)
```

## Arguments

- x:

  A gene-by-cell count matrix, a `SingleCellExperiment`, or a loaded
  Seurat object from an `.rds` file.

- counts_assay:

  Assay name to use when `x` is a `SingleCellExperiment`.

- cell_metadata:

  Optional cell-level metadata with one row per cell.

- feature_metadata:

  Optional feature-level metadata with one row per gene.

- sample_col:

  Optional sample identifier column in `cell_metadata`.

- batch_col:

  Optional batch identifier column in `cell_metadata`.

- condition_col:

  Optional condition column in `cell_metadata`.

- time_col:

  Optional time or ordering column in `cell_metadata`.

- organism:

  Organism label stored in package metadata.

- seurat_assay:

  Assay name to extract when `x` is a Seurat object. If `NULL`, the
  active/default Seurat assay is used.

- seurat_layer:

  Layer to use as raw counts when `x` is a Seurat object. Defaults to
  `"counts"`.

- transfer_reductions:

  Logical indicating whether existing Seurat dimensional reductions (for
  example `pca` and `umap`) should be copied into `reducedDims(sce)`.

## Value

A standardized `SingleCellExperiment`.
