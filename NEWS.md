# scProkaR 0.99.0

## New

* First submission to Bioconductor.
* `CreateBacObject()` builds a `SingleCellExperiment` from a count matrix, a
  10x Genomics directory, or a Seurat object, and `MergeBacObjects()` combines
  several such objects while keeping cell barcodes unique.
* `RunBacQC()` and `FilterBacCells()` compute and apply prokaryote-specific
  quality-control metrics, including ribosomal RNA fraction.
* `IntegrateBacData()`, `RegisterIntegrationEmbedding()`,
  `RunIntegratedUMAP()` and `RunIntegratedClustering()` provide batch
  integration and downstream embedding, with `BenchmarkIntegration()` scoring
  the corrected embeddings using scIB-inspired metrics.
* `aggregate_pseudobulk()`, `normalize_pseudobulk()` and the `run_edger_*()`
  functions provide pseudobulk differential expression built on edgeR,
  covering both pairwise contrasts and spline-based time courses.
* `run_tata()` implements Time Aware Trajectory Analysis (TATA), a
  time-informed trajectory abstraction for perturbation time courses, together
  with branch probability, gene trend, plotting and benchmarking helpers.
* The packaged `sce1` dataset provides a small bacterial single-cell example
  used throughout the vignettes.
