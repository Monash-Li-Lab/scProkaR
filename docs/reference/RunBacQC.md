# Compute bacterial QC metrics

Adds Seurat-style per-cell feature fractions to a
`SingleCellExperiment`, including rRNA and ribosomal-protein content.

## Usage

``` r
RunBacQC(
  sce,
  rrna_pattern = "^(rrs|rrl|rrf)",
  ribo_pattern = "^(rpl|rps)",
  gene_class_col = NULL,
  store = TRUE
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- rrna_pattern:

  Regular expression used to identify rRNA features.

- ribo_pattern:

  Regular expression used to identify ribosomal protein features.

- gene_class_col:

  Optional `rowData` column describing feature classes. When present,
  values matching `"rrna"` are used for rRNA metrics and values matching
  `"ribo"` or `"ribosomal_protein"` are used for ribosomal protein
  metrics.

- store:

  If `TRUE`, store metrics in `colData(sce)` and return the modified
  object. If `FALSE`, return the QC table only.

## Value

A `SingleCellExperiment` when `store = TRUE`, otherwise a data frame of
QC metrics.
