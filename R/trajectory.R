#' Run PAGA-inspired and time/dose trajectory analysis
#'
#' `RunTrajectory()` performs two linked analyses on an existing embedding:
#'
#' 1. A PAGA-inspired cluster abstraction that builds a cell-level kNN graph in
#'    the supplied neighbour space and summarizes cluster-cluster connectivity.
#' 2. A Squidiff-style downstream trajectory summary that works on real cells,
#'    computes per-panel density on the same embedding, tracks time-ordered
#'    centroids, and summarizes cluster composition over time and optionally
#'    dose.
#'
#' The graph abstraction is intentionally PAGA-inspired rather than an exact
#' Scanpy reimplementation. Cluster connectivity is normalized as
#' `observed_edges / min(k * n_i, k * n_j)` and therefore bounded in `[0, 1]`.
#'
#' When both `time_col` and `dose_col` are supplied, dose-specific trajectory
#' panels use the following baseline rule:
#'
#' - for `dose == baseline_dose`, use cells from that dose at all valid times
#' - for `dose > baseline_dose`, use shared baseline cells at
#'   `(baseline_dose, baseline_time)` plus cells from the target dose at
#'   non-baseline timepoints
#'
#' Results are stored in `metadata(object)$SCProkaR$trajectory`, with
#' `paga`, `density`, `composition`, and `params` entries.
#'
#' @param object A `SingleCellExperiment`.
#' @param reduction Candidate reduced-dimension names to use directly as the
#'   neighbour and plotting space. Existing embeddings are preferred and no PCA
#'   is recomputed when the reduction exists.
#' @param cluster_col Metadata column containing existing cluster labels.
#' @param time_col Optional metadata column containing time.
#' @param dose_col Optional metadata column containing dose.
#' @param trajectory_col Optional metadata column used to restrict analysis to
#'   real cells.
#' @param real_value Value in `trajectory_col` indicating real cells.
#' @param n_neighbors Number of neighbours used to build the cell-level graph.
#' @param edge_threshold Minimum connectivity weight to emphasize in plots.
#' @param valid_times Optional explicit ordered time values.
#' @param valid_doses Optional explicit ordered dose values.
#' @param baseline_dose Baseline dose value used for shared-baseline panels.
#' @param baseline_time Baseline time value used for shared-baseline panels.
#' @param density_floor Lower clamp used by density plots.
#' @param kde_n Grid size passed to `MASS::kde2d()` when available.
#' @param min_cells_density Minimum panel size before using a constant fallback
#'   density.
#' @param seed Random seed used for reproducible helper operations.
#' @param verbose Whether to emit warnings for skipped panels and fallbacks.
#' @param dims Optional dimensions of the embedding to use. When `NULL`, all
#'   available dimensions are used for neighbours and the first two are used for
#'   plotting.
#'
#' @return A `SingleCellExperiment` with trajectory results stored in metadata.
#' @export
RunTrajectory <- function(
    object,
    reduction = c("fdl", "X_fdl", "umap", "UMAP", "pca", "PCA"),
    cluster_col = "seurat_clusters",
    time_col = NULL,
    dose_col = NULL,
    trajectory_col = NULL,
    real_value = "real",
    n_neighbors = 15,
    edge_threshold = 0.1,
    valid_times = NULL,
    valid_doses = NULL,
    baseline_dose = 0,
    baseline_time = 0,
    density_floor = 0.2,
    kde_n = 100,
    min_cells_density = 10,
    seed = 1,
    verbose = TRUE,
    dims = NULL
) {
  set.seed(seed)

  meta <- .scprokar_get_cell_metadata(object)
  if (!cluster_col %in% colnames(meta)) {
    stop("Missing cluster column: ", cluster_col, call. = FALSE)
  }

  reduction_name <- .scprokar_resolve_reduction_name(object, reduction)
  emb <- .scprokar_get_embedding_matrix(object, reduction_name, dims = dims)
  xy <- emb[, seq_len(2), drop = FALSE]
  colnames(xy) <- c("x", "y")

  clusters <- .scprokar_cluster_factor(meta[[cluster_col]])
  knn_edges <- .scprokar_compute_knn_edges(emb, k = n_neighbors)
  paga_edges <- .scprokar_aggregate_cluster_connectivity(knn_edges, clusters, k = n_neighbors)
  paga_nodes <- .scprokar_compute_cluster_centroids(xy, clusters)

  density_results <- NULL
  composition <- NULL
  real_mask <- rep(TRUE, nrow(meta))

  if (!is.null(trajectory_col) && trajectory_col %in% colnames(meta)) {
    real_mask <- as.character(meta[[trajectory_col]]) == as.character(real_value)
  }

  if (!is.null(time_col)) {
    if (!time_col %in% colnames(meta)) {
      stop("Missing time column: ", time_col, call. = FALSE)
    }

    use_dose <- !is.null(dose_col) && dose_col %in% colnames(meta)
    if (!is.null(dose_col) && !dose_col %in% colnames(meta)) {
      if (isTRUE(verbose)) {
        warning("`dose_col` was not found; falling back to time-only trajectory summarization.", call. = FALSE)
      }
      use_dose <- FALSE
    }

    real_meta <- meta[real_mask, , drop = FALSE]
    real_xy <- xy[real_mask, , drop = FALSE]
    real_clusters <- clusters[real_mask]

    valid_times <- .scprokar_resolve_values(real_meta[[time_col]], valid_times, label = "time")
    valid_doses <- if (use_dose) {
      .scprokar_resolve_values(real_meta[[dose_col]], valid_doses, label = "dose")
    } else {
      NULL
    }

    density_results <- .scprokar_density_summary(
      xy = real_xy,
      meta = real_meta,
      clusters = real_clusters,
      time_col = time_col,
      dose_col = if (use_dose) dose_col else NULL,
      valid_times = valid_times,
      valid_doses = valid_doses,
      baseline_dose = baseline_dose,
      baseline_time = baseline_time,
      kde_n = kde_n,
      min_cells_density = min_cells_density,
      density_floor = density_floor,
      verbose = verbose
    )
    density_results$params$reduction <- reduction_name
    density_results$params$trajectory_col <- if (!is.null(trajectory_col) && trajectory_col %in% colnames(meta)) trajectory_col else NULL
    density_results$params$real_value <- real_value

    composition <- .scprokar_composition_summary(
      meta = real_meta,
      clusters = real_clusters,
      time_col = time_col,
      dose_col = if (use_dose) dose_col else NULL,
      valid_times = valid_times,
      valid_doses = valid_doses
    )
  }

  result <- list(
    paga = list(
      nodes = paga_nodes,
      edges = paga_edges,
      params = list(
        reduction = reduction_name,
        dims = if (is.null(dims)) seq_len(ncol(emb)) else dims,
        cluster_col = cluster_col,
        n_neighbors = min(n_neighbors, max(1, nrow(emb) - 1)),
        edge_threshold = edge_threshold
      )
    ),
    density = density_results,
    composition = composition,
    params = list(
      reduction = reduction_name,
      cluster_col = cluster_col,
      time_col = time_col,
      dose_col = if (!is.null(dose_col) && dose_col %in% colnames(meta)) dose_col else NULL,
      trajectory_col = if (!is.null(trajectory_col) && trajectory_col %in% colnames(meta)) trajectory_col else NULL,
      real_value = real_value,
      valid_times = valid_times,
      valid_doses = valid_doses,
      baseline_dose = baseline_dose,
      baseline_time = baseline_time,
      density_floor = density_floor,
      edge_threshold = edge_threshold
    )
  )

  sc_meta <- .scprokar_get_metadata(object)
  sc_meta$trajectory <- result
  if (is.null(sc_meta$paga)) {
    sc_meta$paga <- list()
  }
  sc_meta$paga$paga <- result$paga

  .scprokar_set_metadata(object, sc_meta)
}

#' Backward-compatible PAGA-style wrapper
#'
#' `RunBacPAGA()` is retained for compatibility and now delegates to
#' [RunTrajectory()] while preserving the existing function name.
#'
#' @inheritParams RunTrajectory
#' @param sce A `SingleCellExperiment`.
#' @param root_cluster Retained for backward compatibility and stored in
#'   metadata, but not used directly by the embedding-driven abstraction.
#' @param pseudotime_method Retained for backward compatibility and stored in
#'   metadata, but not used directly by the embedding-driven abstraction.
#' @param graph_name Name used when mirroring the PAGA result under the legacy
#'   metadata path `metadata(sce)$SCProkaR$paga[[graph_name]]`.
#'
#' @return A `SingleCellExperiment` with trajectory results stored in metadata.
#' @export
RunBacPAGA <- function(
    sce,
    cluster_col,
    reduction = c("fdl", "X_fdl", "umap", "UMAP", "pca", "PCA"),
    dims = NULL,
    root_cluster = NULL,
    pseudotime_method = c("paga", "mst"),
    graph_name = "paga",
    time_col = NULL,
    dose_col = NULL,
    trajectory_col = NULL,
    real_value = "real",
    n_neighbors = 15,
    edge_threshold = 0.1,
    valid_times = NULL,
    valid_doses = NULL,
    baseline_dose = 0,
    baseline_time = 0,
    density_floor = 0.2,
    kde_n = 100,
    min_cells_density = 10,
    seed = 1,
    verbose = TRUE
) {
  pseudotime_method <- match.arg(pseudotime_method)
  out <- RunTrajectory(
    object = sce,
    reduction = reduction,
    cluster_col = cluster_col,
    time_col = time_col,
    dose_col = dose_col,
    trajectory_col = trajectory_col,
    real_value = real_value,
    n_neighbors = n_neighbors,
    edge_threshold = edge_threshold,
    valid_times = valid_times,
    valid_doses = valid_doses,
    baseline_dose = baseline_dose,
    baseline_time = baseline_time,
    density_floor = density_floor,
    kde_n = kde_n,
    min_cells_density = min_cells_density,
    seed = seed,
    verbose = verbose,
    dims = dims
  )

  meta <- .scprokar_get_metadata(out)
  meta$trajectory$paga$params$legacy_root_cluster <- root_cluster
  meta$trajectory$paga$params$legacy_pseudotime_method <- pseudotime_method
  if (is.null(meta$paga)) {
    meta$paga <- list()
  }
  meta$paga[[graph_name]] <- meta$trajectory$paga

  .scprokar_set_metadata(out, meta)
}

#' Plot the stored PAGA-inspired abstracted graph
#'
#' Overlays the cluster-cluster abstracted graph on the embedding used by
#' [RunTrajectory()] or [RunBacPAGA()].
#'
#' @param object A `SingleCellExperiment` with stored trajectory results.
#' @param reduction Optional reduced-dimension name to use for the background
#'   cells. Defaults to the stored reduction.
#' @param cluster_col Optional cluster column. Defaults to the stored column.
#' @param edge_threshold Minimum edge weight to draw. Defaults to the stored
#'   threshold.
#' @param point_size Background point size.
#' @param point_alpha Background point alpha.
#' @param node_size Centroid node size.
#' @param label_nodes Whether to draw cluster labels.
#' @param palette Optional named or unnamed vector of cluster colours.
#'
#' @return A `ggplot2` object.
#' @export
PlotTrajectoryPAGA <- function(
    object,
    reduction = NULL,
    cluster_col = NULL,
    edge_threshold = NULL,
    point_size = 0.4,
    point_alpha = 0.7,
    node_size = 3,
    label_nodes = TRUE,
    palette = NULL
) {
  traj <- .scprokar_get_trajectory_results(object)
  paga <- traj$paga
  reduction_name <- if (is.null(reduction)) paga$params$reduction else .scprokar_resolve_reduction_name(object, reduction)
  cluster_col <- if (is.null(cluster_col)) paga$params$cluster_col else cluster_col
  edge_threshold <- if (is.null(edge_threshold)) paga$params$edge_threshold else edge_threshold

  meta <- .scprokar_get_cell_metadata(object)
  if (!cluster_col %in% colnames(meta)) {
    stop("Missing cluster column: ", cluster_col, call. = FALSE)
  }

  xy <- .scprokar_get_embedding_matrix(object, reduction_name)
  xy <- xy[, seq_len(2), drop = FALSE]
  colnames(xy) <- c("x", "y")
  clusters <- .scprokar_cluster_factor(meta[[cluster_col]])

  point_df <- data.frame(
    x = xy[, 1],
    y = xy[, 2],
    cluster = clusters,
    stringsAsFactors = FALSE
  )
  nodes <- paga$nodes
  edges <- paga$edges[paga$edges$weight >= edge_threshold, , drop = FALSE]
  edge_df <- .scprokar_join_edge_coords(edges, nodes)
  colors <- .scprokar_discrete_palette(levels(clusters), palette = palette)

  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = point_df,
      mapping = ggplot2::aes(x = x, y = y, color = cluster),
      size = point_size,
      alpha = point_alpha
    )

  if (nrow(edge_df) > 0) {
    plot <- plot +
      ggplot2::geom_segment(
        data = edge_df,
        mapping = ggplot2::aes(
          x = x_from, y = y_from,
          xend = x_to, yend = y_to,
          linewidth = weight
        ),
        inherit.aes = FALSE,
        color = "grey25",
        alpha = 0.85
      ) +
      ggplot2::scale_linewidth(range = c(0.4, 2.5), guide = "none")
  }

  plot <- plot +
    ggplot2::geom_point(
      data = nodes,
      mapping = ggplot2::aes(x = x, y = y, fill = cluster),
      shape = 21,
      color = "black",
      size = node_size,
      stroke = 0.5,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = NULL, color = "Cluster", fill = "Cluster") +
    ggplot2::theme_minimal()

  if (isTRUE(label_nodes)) {
    plot <- plot +
      ggplot2::geom_text(
        data = nodes,
        mapping = ggplot2::aes(x = x, y = y, label = cluster),
        inherit.aes = FALSE,
        size = 3
      )
  }

  plot
}

#' Plot density-based time/dose trajectories
#'
#' Draws real cells in the background, overlays panel-specific point densities on
#' the same embedding, and connects time-ordered centroids with arrows.
#'
#' @param object A `SingleCellExperiment` with stored trajectory density results.
#' @param reduction Optional reduced-dimension name to use for the background
#'   cells. Defaults to the stored reduction.
#' @param point_size Foreground point size.
#' @param bg_point_size Background point size.
#' @param bg_alpha Background point alpha.
#' @param arrow_size Arrow line width.
#' @param label_size Time label text size.
#' @param palette Optional continuous colour palette.
#'
#' @return A `ggplot2` object.
#' @export
PlotTrajectoryDensity <- function(
    object,
    reduction = NULL,
    point_size = 0.8,
    bg_point_size = 0.25,
    bg_alpha = 0.23,
    arrow_size = 0.5,
    label_size = 3,
    palette = NULL
) {
  traj <- .scprokar_get_trajectory_results(object)
  density <- traj$density
  if (is.null(density) || is.null(density$points) || nrow(density$points) == 0) {
    stop("No trajectory density results are stored. Run RunTrajectory() with `time_col` first.", call. = FALSE)
  }

  params <- density$params
  reduction_name <- if (is.null(reduction)) params$reduction else .scprokar_resolve_reduction_name(object, reduction)
  meta <- .scprokar_get_cell_metadata(object)
  xy <- .scprokar_get_embedding_matrix(object, reduction_name)
  xy <- xy[, seq_len(2), drop = FALSE]
  colnames(xy) <- c("x", "y")

  real_mask <- rep(TRUE, nrow(meta))
  if (!is.null(params$trajectory_col) && params$trajectory_col %in% colnames(meta)) {
    real_mask <- as.character(meta[[params$trajectory_col]]) == as.character(params$real_value)
  }

  bg_df <- data.frame(
    x = xy[real_mask, 1],
    y = xy[real_mask, 2],
    stringsAsFactors = FALSE
  )
  point_df <- density$points
  point_df$density_plot <- pmax(params$density_floor, pmin(1, point_df$density))
  centroid_df <- density$centroids
  centroid_df$time_label <- .scprokar_format_time_label(centroid_df$time)
  arrow_df <- .scprokar_centroid_segments(centroid_df, params$valid_times)

  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = bg_df,
      mapping = ggplot2::aes(x = x, y = y),
      color = "grey80",
      size = bg_point_size,
      alpha = bg_alpha
    ) +
    ggplot2::geom_point(
      data = point_df,
      mapping = ggplot2::aes(x = x, y = y, color = density_plot),
      size = point_size
    )

  if (nrow(arrow_df) > 0) {
    plot <- plot +
      ggplot2::geom_segment(
        data = arrow_df,
        mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        linewidth = arrow_size,
        arrow = grid::arrow(length = grid::unit(0.12, "inches")),
        color = "black"
      )
  }

  plot <- plot +
    ggplot2::geom_point(
      data = centroid_df,
      mapping = ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      size = 2.5,
      shape = 21,
      fill = "white",
      color = "black"
    ) +
    ggplot2::geom_text(
      data = centroid_df,
      mapping = ggplot2::aes(x = x, y = y, label = time_label),
      inherit.aes = FALSE,
      size = label_size,
      nudge_y = 0.05
    ) +
    ggplot2::scale_color_gradientn(
      colours = .scprokar_density_palette(palette),
      limits = c(params$density_floor, 1),
      name = "Density"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal()

  if (length(unique(point_df$panel_dose)) > 1) {
    plot <- plot + ggplot2::facet_wrap(~panel_dose, nrow = 1)
  }

  plot
}

#' Plot cluster composition over time and dose
#'
#' Uses stored trajectory composition summaries to show within-panel cluster
#' fractions over time.
#'
#' @param object A `SingleCellExperiment` with stored trajectory composition.
#' @param cluster_col Optional cluster column, retained for API symmetry.
#' @param time_col Optional time column, retained for API symmetry.
#' @param dose_col Optional dose column, retained for API symmetry.
#' @param palette Optional cluster palette.
#'
#' @return A `ggplot2` object.
#' @export
PlotTrajectoryComposition <- function(
    object,
    cluster_col = NULL,
    time_col = NULL,
    dose_col = NULL,
    palette = NULL
) {
  traj <- .scprokar_get_trajectory_results(object)
  composition <- traj$composition
  if (is.null(composition) || nrow(composition) == 0) {
    stop("No trajectory composition summary is stored. Run RunTrajectory() with `time_col` first.", call. = FALSE)
  }

  params <- traj$params
  cluster_levels <- unique(as.character(composition$cluster))
  colors <- .scprokar_discrete_palette(cluster_levels, palette = palette)
  composition$time <- factor(as.character(composition$time), levels = as.character(params$valid_times))

  plot <- ggplot2::ggplot(
    composition,
    ggplot2::aes(x = time, y = prop, color = cluster, group = cluster)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Time", y = "Cluster fraction", color = "Cluster") +
    ggplot2::theme_minimal()

  if (!all(is.na(composition$dose)) && length(unique(composition$dose)) > 1) {
    plot <- plot + ggplot2::facet_wrap(~dose, nrow = 1)
  }

  plot
}

#' @keywords internal
.scprokar_get_cell_metadata <- function(object) {
  as.data.frame(SummarizedExperiment::colData(object))
}

#' @keywords internal
.scprokar_get_trajectory_results <- function(object) {
  meta <- .scprokar_get_metadata(object)
  if (is.null(meta$trajectory)) {
    stop("No trajectory results are stored in metadata(object)$SCProkaR$trajectory.", call. = FALSE)
  }
  meta$trajectory
}

#' @keywords internal
.scprokar_resolve_reduction_name <- function(object, reduction) {
  available <- SingleCellExperiment::reducedDimNames(object)
  candidates <- reduction[!is.na(reduction)]
  found <- candidates[candidates %in% available]
  if (length(found) == 0) {
    stop(
      "None of the requested reductions were found. Requested: ",
      paste(candidates, collapse = ", "),
      ". Available: ",
      paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  found[[1]]
}

#' @keywords internal
.scprokar_get_embedding_matrix <- function(object, reduction, dims = NULL) {
  if (!reduction %in% SingleCellExperiment::reducedDimNames(object)) {
    stop("Reduction not found: ", reduction, call. = FALSE)
  }
  emb <- as.matrix(SingleCellExperiment::reducedDim(object, reduction))
  if (!is.null(dims)) {
    if (max(dims) > ncol(emb)) {
      stop("Requested dimensions exceed the available columns in reduction '", reduction, "'.", call. = FALSE)
    }
    emb <- emb[, dims, drop = FALSE]
  }
  if (ncol(emb) < 2) {
    stop("Reduction '", reduction, "' must have at least two dimensions for trajectory plotting.", call. = FALSE)
  }
  emb
}

#' @keywords internal
.scprokar_cluster_factor <- function(x) {
  if (is.factor(x)) {
    return(factor(as.character(x), levels = levels(x)))
  }
  levels_x <- unique(as.character(x[!is.na(x)]))
  factor(as.character(x), levels = levels_x)
}

#' @keywords internal
.scprokar_resolve_values <- function(x, valid, label) {
  present <- unique(x[!is.na(x)])
  present <- .scprokar_sort_values(present)

  if (is.null(valid)) {
    return(present)
  }

  keep <- valid[valid %in% present]
  if (length(keep) == 0) {
    stop("No requested ", label, " values are present in the metadata.", call. = FALSE)
  }
  keep
}

#' @keywords internal
.scprokar_sort_values <- function(x) {
  if (is.numeric(x)) {
    sort(unique(as.numeric(x)))
  } else {
    sort(unique(as.character(x)))
  }
}

#' @keywords internal
.scprokar_compute_knn_edges <- function(emb, k = 15) {
  idx <- .scprokar_neighbor_index(emb, k = k)
  if (ncol(idx) == 0) {
    return(data.frame(cell_from = integer(), cell_to = integer()))
  }

  edges <- data.frame(
    cell_from = rep(seq_len(nrow(idx)), each = ncol(idx)),
    cell_to = as.vector(t(idx)),
    stringsAsFactors = FALSE
  )
  edges <- edges[edges$cell_from != edges$cell_to, , drop = FALSE]
  edges$key1 <- pmin(edges$cell_from, edges$cell_to)
  edges$key2 <- pmax(edges$cell_from, edges$cell_to)
  edges <- unique(edges[, c("key1", "key2"), drop = FALSE])
  names(edges) <- c("cell_from", "cell_to")
  edges
}

#' @keywords internal
.scprokar_aggregate_cluster_connectivity <- function(knn_edges, clusters, k) {
  levels_cluster <- levels(clusters)
  cluster_size <- table(clusters)

  if (nrow(knn_edges) == 0) {
    return(data.frame(
      from = character(),
      to = character(),
      observed_edges = numeric(),
      weight = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  edge_clusters <- data.frame(
    from = as.character(clusters[knn_edges$cell_from]),
    to = as.character(clusters[knn_edges$cell_to]),
    stringsAsFactors = FALSE
  )
  edge_clusters <- edge_clusters[edge_clusters$from != edge_clusters$to, , drop = FALSE]
  if (nrow(edge_clusters) == 0) {
    return(data.frame(
      from = character(),
      to = character(),
      observed_edges = numeric(),
      weight = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  edge_clusters$idx_from <- match(edge_clusters$from, levels_cluster)
  edge_clusters$idx_to <- match(edge_clusters$to, levels_cluster)
  swap <- edge_clusters$idx_from > edge_clusters$idx_to
  if (any(swap)) {
    tmp <- edge_clusters$from[swap]
    edge_clusters$from[swap] <- edge_clusters$to[swap]
    edge_clusters$to[swap] <- tmp
  }

  observed <- stats::aggregate(
    rep(1, nrow(edge_clusters)),
    by = list(from = edge_clusters$from, to = edge_clusters$to),
    FUN = sum
  )
  names(observed)[3] <- "observed_edges"

  observed$weight <- vapply(seq_len(nrow(observed)), function(i) {
    n_i <- as.numeric(cluster_size[[observed$from[i]]])
    n_j <- as.numeric(cluster_size[[observed$to[i]]])
    max_possible <- min(k * n_i, k * n_j)
    min(1, observed$observed_edges[i] / max(1, max_possible))
  }, numeric(1))

  observed
}

#' @keywords internal
.scprokar_compute_cluster_centroids <- function(xy, clusters) {
  levels_cluster <- levels(clusters)
  centroids <- lapply(levels_cluster, function(cluster) {
    mask <- clusters == cluster
    data.frame(
      cluster = cluster,
      x = mean(xy[mask, 1]),
      y = mean(xy[mask, 2]),
      n_cells = sum(mask),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, centroids)
}

#' @keywords internal
.scprokar_build_dose_panel_mask <- function(
    meta,
    dose_col,
    time_col,
    target_dose,
    valid_times,
    baseline_dose,
    baseline_time
) {
  dose_values <- meta[[dose_col]]
  time_values <- meta[[time_col]]

  if (identical(as.character(target_dose), as.character(baseline_dose))) {
    mask <- as.character(dose_values) == as.character(baseline_dose) &
      time_values %in% valid_times
    baseline_flag <- as.character(dose_values) == as.character(baseline_dose) &
      as.character(time_values) == as.character(baseline_time)
  } else {
    mask <- (
      as.character(dose_values) == as.character(baseline_dose) &
        as.character(time_values) == as.character(baseline_time)
    ) | (
      as.character(dose_values) == as.character(target_dose) &
        time_values %in% setdiff(valid_times, baseline_time)
    )
    baseline_flag <- as.character(dose_values) == as.character(baseline_dose) &
      as.character(time_values) == as.character(baseline_time)
  }

  list(mask = mask, is_baseline = baseline_flag)
}

#' @keywords internal
.scprokar_density_summary <- function(
    xy,
    meta,
    clusters,
    time_col,
    dose_col = NULL,
    valid_times,
    valid_doses = NULL,
    baseline_dose,
    baseline_time,
    kde_n,
    min_cells_density,
    density_floor,
    verbose
) {
  panel_points <- list()
  panel_centroids <- list()
  panel_index <- 0

  if (is.null(dose_col)) {
    panel_specs <- list(list(panel_dose = "all", dose = NA))
  } else {
    panel_specs <- lapply(valid_doses, function(dose) {
      list(panel_dose = as.character(dose), dose = dose)
    })
  }

  for (spec in panel_specs) {
    panel_index <- panel_index + 1

    if (is.null(dose_col)) {
      mask <- meta[[time_col]] %in% valid_times
      baseline_flag <- as.character(meta[[time_col]]) == as.character(baseline_time)
    } else {
      panel_mask <- .scprokar_build_dose_panel_mask(
        meta = meta,
        dose_col = dose_col,
        time_col = time_col,
        target_dose = spec$dose,
        valid_times = valid_times,
        baseline_dose = baseline_dose,
        baseline_time = baseline_time
      )
      mask <- panel_mask$mask
      baseline_flag <- panel_mask$is_baseline
    }

    if (!any(mask)) {
      if (isTRUE(verbose)) {
        warning("Skipping trajectory panel '", spec$panel_dose, "' because it contains no cells.", call. = FALSE)
      }
      next
    }

    xy_panel <- xy[mask, , drop = FALSE]
    meta_panel <- meta[mask, , drop = FALSE]
    density <- .scprokar_compute_point_density(xy_panel, n = kde_n, min_cells_density = min_cells_density)

    panel_points[[panel_index]] <- data.frame(
      panel_dose = spec$panel_dose,
      x = xy_panel[, 1],
      y = xy_panel[, 2],
      density = density,
      time = meta_panel[[time_col]],
      dose = if (is.null(dose_col)) NA else meta_panel[[dose_col]],
      is_baseline = baseline_flag[mask],
      stringsAsFactors = FALSE
    )

    centroid_rows <- lapply(valid_times, function(time_value) {
      time_mask <- as.character(meta_panel[[time_col]]) == as.character(time_value)
      if (!any(time_mask)) {
        if (isTRUE(verbose)) {
          warning(
            "Skipping centroid for panel '", spec$panel_dose,
            "' at time '", time_value, "' because it contains no cells.",
            call. = FALSE
          )
        }
        return(NULL)
      }
      data.frame(
        panel_dose = spec$panel_dose,
        time = time_value,
        x = mean(xy_panel[time_mask, 1]),
        y = mean(xy_panel[time_mask, 2]),
        n_cells = sum(time_mask),
        stringsAsFactors = FALSE
      )
    })

    centroid_rows <- Filter(Negate(is.null), centroid_rows)
    if (length(centroid_rows) > 0) {
      panel_centroids[[panel_index]] <- do.call(rbind, centroid_rows)
    }
  }

  points <- if (length(panel_points) > 0) do.call(rbind, panel_points) else NULL
  centroids <- if (length(panel_centroids) > 0) do.call(rbind, panel_centroids) else NULL

  list(
    points = points,
    centroids = centroids,
    params = list(
      reduction = NULL,
      time_col = time_col,
      dose_col = dose_col,
      valid_times = valid_times,
      valid_doses = valid_doses,
      baseline_dose = baseline_dose,
      baseline_time = baseline_time,
      density_floor = density_floor
    )
  )
}

#' @keywords internal
.scprokar_compute_point_density <- function(xy, n = 100, min_cells_density = 10) {
  xy <- as.matrix(xy)
  if (nrow(xy) < min_cells_density) {
    return(rep(1, nrow(xy)))
  }

  if (requireNamespace("MASS", quietly = TRUE) &&
      length(unique(xy[, 1])) > 1 &&
      length(unique(xy[, 2])) > 1) {
    kde <- tryCatch(
      MASS::kde2d(xy[, 1], xy[, 2], n = n),
      error = function(e) NULL
    )
    if (!is.null(kde)) {
      raw_density <- .scprokar_interp_kde2d(kde, xy)
      return(.scprokar_normalize_density(raw_density))
    }
  }

  rep(1, nrow(xy))
}

#' @keywords internal
.scprokar_interp_kde2d <- function(kde, xy) {
  ix <- findInterval(xy[, 1], kde$x, all.inside = TRUE)
  iy <- findInterval(xy[, 2], kde$y, all.inside = TRUE)
  vapply(seq_len(nrow(xy)), function(i) kde$z[ix[i], iy[i]], numeric(1))
}

#' @keywords internal
.scprokar_normalize_density <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (isTRUE(all.equal(rng[1], rng[2]))) {
    return(rep(1, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

#' @keywords internal
.scprokar_composition_summary <- function(meta, clusters, time_col, dose_col = NULL, valid_times, valid_doses = NULL) {
  dose_values <- if (is.null(dose_col)) rep("all", nrow(meta)) else meta[[dose_col]]
  out <- expand.grid(
    dose = if (is.null(dose_col)) "all" else as.character(valid_doses),
    time = as.character(valid_times),
    cluster = levels(clusters),
    stringsAsFactors = FALSE
  )

  out$n_cells <- mapply(function(dose, time, cluster) {
    mask <- as.character(meta[[time_col]]) == as.character(time) &
      as.character(clusters) == as.character(cluster)
    if (!is.null(dose_col)) {
      mask <- mask & as.character(dose_values) == as.character(dose)
    }
    sum(mask)
  }, out$dose, out$time, out$cluster)

  totals <- stats::aggregate(n_cells ~ dose + time, data = out, FUN = sum)
  names(totals)[3] <- "panel_total"
  out <- merge(out, totals, by = c("dose", "time"), all.x = TRUE, sort = FALSE)
  out$prop <- ifelse(out$panel_total > 0, out$n_cells / out$panel_total, 0)
  out <- out[, c("dose", "time", "cluster", "n_cells", "prop")]
  out
}

#' @keywords internal
.scprokar_join_edge_coords <- function(edges, nodes) {
  if (nrow(edges) == 0) {
    return(data.frame(
      from = character(), to = character(), observed_edges = numeric(), weight = numeric(),
      x_from = numeric(), y_from = numeric(), x_to = numeric(), y_to = numeric()
    ))
  }

  node_from <- nodes[, c("cluster", "x", "y")]
  names(node_from) <- c("from", "x_from", "y_from")
  node_to <- nodes[, c("cluster", "x", "y")]
  names(node_to) <- c("to", "x_to", "y_to")
  merge(merge(edges, node_from, by = "from", sort = FALSE), node_to, by = "to", sort = FALSE)
}

#' @keywords internal
.scprokar_centroid_segments <- function(centroids, valid_times) {
  if (is.null(centroids) || nrow(centroids) < 2) {
    return(data.frame())
  }

  centroids$time_factor <- match(as.character(centroids$time), as.character(valid_times))
  split_panels <- split(centroids, centroids$panel_dose)
  segments <- lapply(split_panels, function(df) {
    df <- df[order(df$time_factor), , drop = FALSE]
    if (nrow(df) < 2) {
      return(NULL)
    }
    data.frame(
      panel_dose = df$panel_dose[-nrow(df)],
      x = df$x[-nrow(df)],
      y = df$y[-nrow(df)],
      xend = df$x[-1],
      yend = df$y[-1],
      stringsAsFactors = FALSE
    )
  })
  segments <- Filter(Negate(is.null), segments)
  if (length(segments) == 0) {
    return(data.frame())
  }
  do.call(rbind, segments)
}

#' @keywords internal
.scprokar_discrete_palette <- function(levels_vec, palette = NULL) {
  if (!is.null(palette)) {
    if (is.null(names(palette))) {
      names(palette) <- levels_vec[seq_len(min(length(levels_vec), length(palette)))]
    }
    return(palette[levels_vec])
  }
  colors <- grDevices::hcl.colors(length(levels_vec), "Dark 3")
  stats::setNames(colors, levels_vec)
}

#' @keywords internal
.scprokar_density_palette <- function(palette = NULL) {
  if (!is.null(palette)) {
    return(palette)
  }
  c("#d9ecff", "#90d4c5", "#f2e96b", "#f6a04d", "#cc4c3c")
}

#' @keywords internal
.scprokar_format_time_label <- function(x) {
  if (is.numeric(x)) {
    paste0(x, "h")
  } else {
    as.character(x)
  }
}
