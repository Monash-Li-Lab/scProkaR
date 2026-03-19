#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SCProkaR)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(ggplot2)
})

# SCProkaR two-sample workflow
# ----------------------------
# Edit the paths and metadata section below, then run:
# Rscript inst/scripts/scprokar_integration_workflow.R

sample1_rds <- "sample1.rds"
sample2_rds <- "sample2.rds"

out_dir <- "scprokar_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Loading Seurat .rds files")
seu1 <- readRDS(sample1_rds)
seu2 <- readRDS(sample2_rds)

# Add sample origin columns if they are not already present.
# Replace these with your own metadata columns if needed.
seu1$sample_id <- "sample1"
seu1$batch <- "sample1"

seu2$sample_id <- "sample2"
seu2$batch <- "sample2"

message("Creating SCProkaR objects")
sce1 <- CreateBacObject(
  seu1,
  seurat_assay = "RNA",
  seurat_layer = "counts",
  sample_col = "sample_id",
  batch_col = "batch"
)

sce2 <- CreateBacObject(
  seu2,
  seurat_assay = "RNA",
  seurat_layer = "counts",
  sample_col = "sample_id",
  batch_col = "batch"
)

message("Running QC")
sce1 <- RunBacQC(sce1)
sce2 <- RunBacQC(sce2)

write.csv(
  as.data.frame(SummarizedExperiment::colData(sce1)),
  file.path(out_dir, "sample1_qc_metrics.csv"),
  row.names = TRUE
)
write.csv(
  as.data.frame(SummarizedExperiment::colData(sce2)),
  file.path(out_dir, "sample2_qc_metrics.csv"),
  row.names = TRUE
)

message("Filtering cells")
sce1 <- FilterBacCells(
  sce1,
  min_counts = 200,
  min_features = 50,
  max_rrna_fraction = 0.20
)

sce2 <- FilterBacCells(
  sce2,
  min_counts = 200,
  min_features = 50,
  max_rrna_fraction = 0.20
)

message("Merging samples")
sce <- MergeBacObjects(sce1, sce2, gene_mode = "intersect")

cluster_col <- if ("seurat_clusters" %in% colnames(SummarizedExperiment::colData(sce))) {
  "seurat_clusters"
} else if ("cluster" %in% colnames(SummarizedExperiment::colData(sce))) {
  "cluster"
} else {
  NULL
}

message("Plotting unintegrated reductions")
p_unint_batch <- PlotReduction(sce, reduction = "umap", colour_by = "batch")
ggsave(
  filename = file.path(out_dir, "01_unintegrated_umap_by_batch.pdf"),
  plot = p_unint_batch,
  width = 6,
  height = 5
)

if (!is.null(cluster_col)) {
  p_unint_cluster <- PlotReduction(sce, reduction = "umap", colour_by = cluster_col)
  ggsave(
    filename = file.path(out_dir, "02_unintegrated_umap_by_cluster.pdf"),
    plot = p_unint_cluster,
    width = 6,
    height = 5
  )
}

message("Running batch integration")
# Choose one method:
integration_method <- "mnn"
# integration_method <- "harmony"

sce <- IntegrateBacData(
  sce,
  batch_col = "batch",
  method = integration_method,
  dims = 1:30
)

integrated_name <- paste0("integrated_", integration_method)

p_integrated_latent_batch <- PlotReduction(sce, reduction = integrated_name, colour_by = "batch")
ggsave(
  filename = file.path(out_dir, "03_integrated_latent_by_batch.pdf"),
  plot = p_integrated_latent_batch,
  width = 6,
  height = 5
)

if (!is.null(cluster_col)) {
  p_integrated_latent_cluster <- PlotReduction(sce, reduction = integrated_name, colour_by = cluster_col)
  ggsave(
    filename = file.path(out_dir, "04_integrated_latent_by_cluster.pdf"),
    plot = p_integrated_latent_cluster,
    width = 6,
    height = 5
  )
}

message("Running UMAP on integrated embedding")
sce <- RunIntegratedUMAP(sce, reduction = integrated_name)
integrated_umap_name <- paste0("umap_", integrated_name)

p_integrated_umap_batch <- PlotReduction(sce, reduction = integrated_umap_name, colour_by = "batch")
ggsave(
  filename = file.path(out_dir, "05_integrated_umap_by_batch.pdf"),
  plot = p_integrated_umap_batch,
  width = 6,
  height = 5
)

if (!is.null(cluster_col)) {
  p_integrated_umap_cluster <- PlotReduction(sce, reduction = integrated_umap_name, colour_by = cluster_col)
  ggsave(
    filename = file.path(out_dir, "06_integrated_umap_by_cluster.pdf"),
    plot = p_integrated_umap_cluster,
    width = 6,
    height = 5
  )
}

message("Clustering cells on integrated embedding")
sce <- RunIntegratedClustering(
  sce,
  reduction = integrated_name,
  cluster_col = "integrated_clusters",
  k = 20,
  resolution = 0.8
)

p_integrated_clusters <- PlotReduction(
  sce,
  reduction = integrated_umap_name,
  colour_by = "integrated_clusters"
)
ggsave(
  filename = file.path(out_dir, "07_integrated_umap_by_integrated_clusters.pdf"),
  plot = p_integrated_clusters,
  width = 6,
  height = 5
)

p_integrated_clusters_split <- PlotReduction(
  sce,
  reduction = integrated_umap_name,
  colour_by = "integrated_clusters",
  facet_by = "batch"
)
ggsave(
  filename = file.path(out_dir, "08_integrated_umap_by_integrated_clusters_split.pdf"),
  plot = p_integrated_clusters_split,
  width = 10,
  height = 4
)

message("Generating Seurat-style integration summary panels")
overview_plots <- PlotIntegrationOverview(
  sce,
  unintegrated_reduction = "umap",
  integrated_reduction = integrated_umap_name,
  batch_col = "batch",
  label_col = "integrated_clusters",
  split_by = "batch"
)

ggsave(
  filename = file.path(out_dir, "09_overview_unintegrated_by_batch.pdf"),
  plot = overview_plots$unintegrated_by_batch,
  width = 6,
  height = 5
)
ggsave(
  filename = file.path(out_dir, "10_overview_integrated_by_batch.pdf"),
  plot = overview_plots$integrated_by_batch,
  width = 6,
  height = 5
)

if (!is.null(overview_plots$unintegrated_by_label)) {
  ggsave(
    filename = file.path(out_dir, "11_overview_unintegrated_by_label.pdf"),
    plot = overview_plots$unintegrated_by_label,
    width = 6,
    height = 5
  )
}
if (!is.null(overview_plots$integrated_by_label)) {
  ggsave(
    filename = file.path(out_dir, "12_overview_integrated_by_label.pdf"),
    plot = overview_plots$integrated_by_label,
    width = 6,
    height = 5
  )
}
if (!is.null(overview_plots$integrated_split)) {
  ggsave(
    filename = file.path(out_dir, "13_overview_integrated_split.pdf"),
    plot = overview_plots$integrated_split,
    width = 10,
    height = 4
  )
}

if (!is.null(cluster_col)) {
  message("Running integration benchmark")
  bench <- BenchmarkIntegration(
    sce,
    batch_col = "batch",
    label_col = cluster_col,
    methods = integrated_name
  )

  write.csv(
    bench$ranking,
    file.path(out_dir, "integration_benchmark_ranking.csv"),
    row.names = FALSE
  )

  ggsave(
    filename = file.path(out_dir, "14_integration_benchmark_overall.pdf"),
    plot = bench$plots$overall,
    width = 6,
    height = 4
  )

  ggsave(
    filename = file.path(out_dir, "15_integration_benchmark_heatmap.pdf"),
    plot = bench$plots$heatmap,
    width = 7,
    height = 4
  )
}

saveRDS(sce, file = file.path(out_dir, "scprokar_integrated_sce.rds"))
message("Workflow complete. Outputs saved to: ", normalizePath(out_dir))
