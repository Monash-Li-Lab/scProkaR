
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCProkaR

<!-- badges: start -->

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://monash-li-lab.github.io/scProka/)
<!-- badges: end -->

SCProkaR is a bacterial and microbial single-cell RNA-seq toolkit built
around `SingleCellExperiment`. It standardizes raw inputs, calculates
bacterial QC metrics such as rRNA fraction, integrates batches with
lightweight backends, benchmarks corrected embeddings including
user-supplied embeddings, and supports integrated clustering for
downstream analysis. The package also includes pseudobulk differential
expression utilities and Time Aware Trajectory Analysis (TATA) for
time-informed trajectory inference.

## Installation

You can install the development version of SCProkaR like so:

``` r
remotes::install_github("Monash-Li-Lab/scProka")
```

## Package website

- Core workflow vignette:
  <https://monash-li-lab.github.io/scProka/articles/scprokar-workflow.html>
- Comprehensive PDF tutorial source:
  `vignettes/scprokar-pdf-tutorial.Rmd`
- TATA vignette:
  <https://monash-li-lab.github.io/scProka/articles/tata-workflow.html>
- Function reference:
  <https://monash-li-lab.github.io/scProka/reference/index.html>

## Example data and workflow

The package ships with `sce1`, a small three-treatment
`SingleCellExperiment` example object used throughout the main workflow
vignette.

``` r
library(SCProkaR)
data("sce1", package = "SCProkaR")

sce1
#> class: SingleCellExperiment 
#> dim: 3722 5000 
#> metadata(0):
#> assays(1): counts
#> rownames(3722): dnaA ABUW-RS00010 ... ABUW-RS20075 ABUW-RS20460
#> rowData names(0):
#> colnames(5000): P1_1_GGAGCTGAGAGGTCCGTCTA P1_1_TATTGTCATCGGAGCGCTGC ...
#>   control4_TCGCGAAGCAGAGCTGCAAC control4_CGGATGCGAACGATTGATGA
#> colData names(4): sample clusters treatments timepoints
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
table(sce1$treatments)
#> 
#> control  PMB0.5    PMB2 
#>    2000    1500    1500
table(sce1$timepoints)
#> 
#>    0    1    4    7 
#>  500 1500 1500 1500
```

A minimal workflow is:

``` r
sce <- sce1
sce$batch <- factor(sce$treatments)
sce <- RunBacQC(sce)

if (requireNamespace("batchelor", quietly = TRUE)) {
  sce <- IntegrateBacData(sce, batch_col = "batch", method = "mnn", dims = 1:10)
  integrated_reduction <- "integrated_mnn"
} else {
  fallback_latent <- matrix(
    rnorm(ncol(sce) * 10),
    ncol = 10,
    dimnames = list(colnames(sce), paste0("dim", seq_len(10)))
  )
  sce <- RegisterIntegrationEmbedding(sce, embedding = fallback_latent, method_name = "mock")
  integrated_reduction <- "integrated_mock"
}

sce <- RunIntegratedUMAP(sce, reduction = integrated_reduction)
sce <- RunIntegratedClustering(sce, reduction = integrated_reduction, cluster_col = "integrated_clusters")

tata <- run_tata(
  sce = sce,
  dimred = integrated_reduction,
  time_col = "timepoints",
  cluster_col = "tata_cluster",
  k = 20,
  do_pseudotime = TRUE,
  seed = 101
)

names(tata)
#>  [1] "sce"                    "cell_graph"             "adjacency_matrix"      
#>  [4] "cluster_graph"          "cluster_edge_table"     "cluster_vertex_table"  
#>  [7] "cluster_pseudotime"     "cell_pseudotime"        "cell_pseudotime_scaled"
#> [10] "parameters"             "topology_table"         "time_table"
```

If your input is a Seurat `.rds` object, you can pass it directly to
`CreateBacObject()` and SCProkaR will extract the assay counts,
metadata, and stored reductions:

``` r
seu <- readRDS("sample1.rds")
sce <- CreateBacObject(
  seu,
  seurat_assay = "RNA",
  seurat_layer = "counts",
  sample_col = "orig.ident"
)
```

For multiple samples, create and QC each sample separately, then merge
with `MergeBacObjects()` before integration:

``` r
seu1 <- readRDS("sample1.rds")
seu2 <- readRDS("sample2.rds")

seu1$sample_id <- "sample1"
seu1$batch <- "sample1"
seu2$sample_id <- "sample2"
seu2$batch <- "sample2"

sce1 <- CreateBacObject(seu1, seurat_assay = "RNA", sample_col = "sample_id", batch_col = "batch")
sce2 <- CreateBacObject(seu2, seurat_assay = "RNA", sample_col = "sample_id", batch_col = "batch")

sce1 <- RunBacQC(sce1)
sce2 <- RunBacQC(sce2)

sce1 <- FilterBacCells(sce1, min_counts = 200, min_features = 50, max_rrna_fraction = 0.2)
sce2 <- FilterBacCells(sce2, min_counts = 200, min_features = 50, max_rrna_fraction = 0.2)

sce <- MergeBacObjects(sce1, sce2, gene_mode = "intersect")
sce <- IntegrateBacData(sce, batch_col = "batch", method = "mnn")
sce <- RunIntegratedUMAP(sce, reduction = "integrated_mnn")
sce <- RunIntegratedClustering(sce, reduction = "integrated_mnn", cluster_col = "integrated_clusters")

PlotReduction(sce, reduction = "umap", colour_by = "batch")
PlotReduction(sce, reduction = "umap_integrated_mnn", colour_by = "batch")
PlotReduction(sce, reduction = "umap_integrated_mnn", colour_by = "integrated_clusters")
```

If the original cell barcodes overlap across samples,
`MergeBacObjects()` will automatically make them unique by prefixing
them with the sample or batch label when available.

## Optional Backends

Some workflows require optional packages:

- `batchelor` for MNN integration
- `harmony` for Harmony integration
- `SeuratObject` only when you pass a Seurat `.rds` directly into
  `CreateBacObject()`
- `igraph` for integrated graph clustering
- `edgeR` for pseudobulk differential expression
- `pkgdown` for building the package website

The public API is intentionally centered on one object type so
downstream methods can share assays, metadata, reduced dimensions, and
analysis results. Custom embeddings stored in `reducedDims(sce)` can
also be benchmarked directly.
