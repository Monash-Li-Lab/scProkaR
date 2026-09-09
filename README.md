# SCProkaR

SCProkaR is a `SingleCellExperiment`-centred toolkit for microbial single-cell
RNA sequencing. It standardises bacterial input objects, computes
prokaryote-specific quality-control metrics such as ribosomal RNA fraction,
integrates batches with several backends, and scores the corrected embeddings
with scIB-inspired benchmarking metrics. It also provides pseudobulk
differential expression built on edgeR, and Time Aware Trajectory Analysis
(TATA), a time-informed trajectory abstraction for perturbation time courses.

## Installation

SCProkaR is under review for Bioconductor. Once it is accepted, install it
with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SCProkaR")
```

To install the development version from GitHub:

```r
BiocManager::install("Monash-Li-Lab/SCProkaR")
```

## Getting started

```r
library(SCProkaR)
data("sce1", package = "SCProkaR")
sce1
```

The package ships `sce1`, a bacterial single-cell time course of 5000 cells
across an untreated control and two polymyxin B concentrations, sampled at
0, 1, 4 and 7 hours.

A minimal workflow:

```r
sce <- sce1
sce$batch <- factor(sce$treatments)
sce <- RunBacQC(sce)
sce <- FilterBacCells(sce, min_counts = 200, min_features = 50)

sce <- IntegrateBacData(sce, batch_col = "batch", method = "mnn", dims = 1:10)
sce <- RunIntegratedUMAP(sce, reduction = "integrated_mnn")
sce <- RunIntegratedClustering(sce, reduction = "integrated_mnn")

PlotReduction(sce, reduction = "umap_integrated_mnn", colour_by = "batch")
```

Time-aware trajectory analysis on the integrated embedding:

```r
set.seed(101)
tata <- run_tata(
    sce,
    dimred = "integrated_mnn",
    time_col = "timepoints",
    cluster_col = "tata_cluster",
    k = 20
)
plot_tata_cluster_graph(tata)
```

## Vignettes

Two vignettes document the package in full:

```r
browseVignettes("SCProkaR")
```

* **SCProkaR core workflow** covers object creation, quality control,
  integration, benchmarking, clustering and pseudobulk differential
  expression.
* **Time Aware Trajectory Analysis** covers the TATA graph, pseudotime,
  branch probabilities and gene trends.

## Main functions

| Stage | Functions |
| --- | --- |
| Object creation | `CreateBacObject()`, `MergeBacObjects()` |
| Quality control | `RunBacQC()`, `FilterBacCells()` |
| Integration | `IntegrateBacData()`, `RegisterIntegrationEmbedding()`, `RunIntegratedUMAP()`, `RunIntegratedClustering()` |
| Benchmarking | `BenchmarkIntegration()` |
| Visualisation | `PlotReduction()`, `PlotIntegrationOverview()` |
| Pseudobulk DE | `aggregate_pseudobulk()`, `normalize_pseudobulk()`, `run_edger_pairwise_de()`, `run_edger_spline_de()`, `run_scproka_de()` |
| Trajectory | `run_tata()`, `compute_tata_pseudotime()`, `compute_tata_branch_probabilities()`, `compute_tata_gene_trends()` |

## Optional backends

Several backends are optional and are only needed for specific steps:

* `batchelor` for mutual nearest neighbour integration
* `harmony` for Harmony integration
* `uwot` for UMAP
* `SeuratObject` for reading Seurat objects into `CreateBacObject()`
* `DropletUtils` for reading 10x Genomics directories
* `scran`, `irlba`, `BiocSingular`, `BiocNeighbors`, `RANN`, `FNN`, `cluster`
  for individual preprocessing and neighbour-search steps

## Getting help

Please open an issue at
<https://github.com/Monash-Li-Lab/SCProkaR/issues>, or ask on the
[Bioconductor support site](https://support.bioconductor.org/) using the
`SCProkaR` tag.

## License

MIT. See `LICENSE`.
