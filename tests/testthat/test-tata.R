test_that("run_tata returns graph, pseudotime, and plotting outputs", {
  sce <- make_toy_multidrug_sce(n_cells = 900, n_features = 70, seed = 31)

  tata <- run_tata(
    sce = sce,
    dimred = "PCA",
    time_col = "timepoint",
    cluster_col = "tata_cluster",
    k = 20,
    cluster_method = "louvain",
    prune_threshold = 0.001,
    time_weight_mode = "directional_confidence",
    seed = 31
  )

  expect_true(all(c(
    "sce",
    "cell_graph",
    "adjacency_matrix",
    "cluster_graph",
    "cluster_edge_table",
    "cluster_pseudotime",
    "cell_pseudotime"
  ) %in% names(tata)))
  expect_s4_class(tata$sce, "SingleCellExperiment")
  expect_true(inherits(tata$cell_graph, "igraph"))
  expect_true(inherits(tata$cluster_graph, "igraph"))
  expect_true(nrow(tata$cluster_edge_table) > 0)
  expect_true("tata_pseudotime" %in% colnames(SummarizedExperiment::colData(tata$sce)))

  graph_plot <- plot_tata_cluster_graph(
    tata,
    dimred = "UMAP",
    layout_mode = "embedding",
    show_cells = TRUE,
    colour_by = "timepoint",
    node_colour_by = "cluster"
  )
  trajectory_plot <- plot_tata_trajectory_embedding(
    tata,
    dimred = "UMAP",
    colour_by = "timepoint",
    mode = "combined",
    max_paths = 3
  )

  expect_s3_class(graph_plot, "ggplot")
  expect_s3_class(trajectory_plot, "ggplot")
})

test_that("TATA benchmark helpers return finite summary statistics", {
  sce <- make_toy_multidrug_sce(n_cells = 900, n_features = 70, seed = 32)
  tata <- run_tata(
    sce = sce,
    dimred = "PCA",
    time_col = "timepoint",
    cluster_col = "tata_cluster",
    k = 20,
    cluster_method = "louvain",
    prune_threshold = 0.001,
    time_weight_mode = "directional_confidence",
    seed = 32
  )

  local_metrics <- local_temporal_order_metrics(
    pseudotime = tata$cell_pseudotime,
    adjacency = tata$adjacency_matrix,
    timepoint = tata$sce$timepoint
  )
  auc_metrics <- adjacent_timepoint_auc(
    pseudotime = tata$cell_pseudotime,
    timepoint = tata$sce$timepoint,
    group = tata$sce$condition,
    min_cells = 5
  )
  onset_metrics <- branch_onset_metrics(
    pseudotime = tata$cell_pseudotime,
    branch = as.character(tata$sce$condition),
    branch_activation = tata$sce$branch_activation,
    true_onset = c(drug_A = 24, drug_C = 48, drug_B = 96)
  )
  reachability <- directed_reachability_concordance(
    tata$cluster_graph,
    clusters = tata$sce$tata_cluster,
    timepoint = tata$sce$timepoint
  )
  causal_mass <- anticausal_edge_mass(
    tata$cluster_graph,
    clusters = tata$sce$tata_cluster,
    timepoint = tata$sce$timepoint
  )

  expect_true(is.data.frame(local_metrics))
  expect_true(is.list(auc_metrics))
  expect_true(is.data.frame(onset_metrics$summary))
  expect_true(is.data.frame(reachability$summary))
  expect_true(is.data.frame(causal_mass$summary))
  expect_true(is.finite(local_metrics$local_temporal_direction_score))
  expect_true(is.finite(auc_metrics$summary$adjacent_timepoint_auc))
})
