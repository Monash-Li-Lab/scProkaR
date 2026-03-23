#' Compute TATA pseudotime
#'
#' Compute a simple cluster-level pseudotime on the directed TATA graph using
#' shortest-path distances from the root cluster with inverse edge weights, then
#' propagate that pseudotime to cells using within-cluster time ranks. If some
#' clusters are not reachable in the directed graph, TATA falls back to the
#' undirected abstract topology for those clusters so late branches are not
#' dropped simply because an intermediate direction is ambiguous.
#'
#' @param cluster_graph Directed cluster-level TATA graph.
#' @param sce A `SingleCellExperiment`.
#' @param cluster_col Column in `colData(sce)` containing inferred clusters.
#' @param time_col Column in `colData(sce)` containing experimental time labels.
#' @param root_cluster Optional root cluster. If `NULL`, the root is chosen as
#'   the cluster with the earliest median timepoint.
#' @param within_cluster_fraction Scale factor used to spread cells within each
#'   cluster after assigning cluster-level pseudotime.
#' @param eps Small constant to avoid division by zero.
#'
#' @return A list with root cluster, cluster pseudotime table, raw cell
#'   pseudotime, scaled cell pseudotime, and the within-cluster step size.
#' @export
compute_tata_pseudotime <- function(
    cluster_graph,
    sce,
    cluster_col = "cluster",
    time_col = "timepoint",
    root_cluster = NULL,
    within_cluster_fraction = 0.5,
    eps = 1e-8) {
  if (!inherits(cluster_graph, "igraph")) {
    stop("`cluster_graph` must be an igraph object.", call. = FALSE)
  }

  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  if (!cluster_col %in% colnames(meta)) {
    stop("Cluster column `", cluster_col, "` was not found in colData.", call. = FALSE)
  }
  if (!time_col %in% colnames(meta)) {
    stop("Time column `", time_col, "` was not found in colData.", call. = FALSE)
  }

  clusters <- factor(meta[[cluster_col]], levels = igraph::V(cluster_graph)$name)
  time_numeric <- .coerce_time_to_numeric(meta[[time_col]])

  cluster_size <- tapply(seq_along(clusters), clusters, length)
  cluster_median_time <- tapply(time_numeric, clusters, stats::median)

  if (is.null(root_cluster)) {
    earliest_time <- min(time_numeric, na.rm = TRUE)
    earliest_fraction <- tapply(time_numeric == earliest_time, clusters, mean)
    earliest_fraction[is.na(earliest_fraction)] <- 0

    candidate_mask <- earliest_fraction > 0
    candidates <- names(earliest_fraction)[candidate_mask]
    if (length(candidates) == 0L) {
      earliest_cluster_time <- min(cluster_median_time, na.rm = TRUE)
      candidates <- names(cluster_median_time)[cluster_median_time == earliest_cluster_time]
    }

    undirected_graph <- igraph::as_undirected(cluster_graph, mode = "collapse")
    edge_weights <- igraph::E(undirected_graph)$weight
    if (length(edge_weights) > 0L) {
      closeness_score <- igraph::closeness(
        undirected_graph,
        mode = "all",
        weights = 1 / pmax(edge_weights, eps),
        normalized = TRUE
      )
    } else {
      closeness_score <- rep(0, igraph::vcount(undirected_graph))
      names(closeness_score) <- igraph::V(undirected_graph)$name
    }

    candidate_table <- data.frame(
      cluster = candidates,
      earliest_fraction = as.numeric(earliest_fraction[candidates]),
      closeness = as.numeric(closeness_score[candidates]),
      size = as.numeric(cluster_size[candidates]),
      stringsAsFactors = FALSE
    )
    candidate_table <- candidate_table[
      order(
        -candidate_table$earliest_fraction,
        -candidate_table$closeness,
        -candidate_table$size
      ),
      ,
      drop = FALSE
    ]

    root_cluster <- candidate_table$cluster[1]
  }

  if (!root_cluster %in% igraph::V(cluster_graph)$name) {
    stop("`root_cluster` is not present in the cluster graph.", call. = FALSE)
  }

  if (igraph::ecount(cluster_graph) == 0L) {
    directed_distance <- rep(Inf, igraph::vcount(cluster_graph))
    names(directed_distance) <- igraph::V(cluster_graph)$name
    directed_distance[root_cluster] <- 0
    undirected_distance <- directed_distance
  } else {
    inverse_weight <- 1 / pmax(igraph::E(cluster_graph)$weight, eps)
    directed_distance <- as.numeric(
      igraph::distances(
        cluster_graph,
        v = root_cluster,
        to = igraph::V(cluster_graph),
        mode = "out",
        weights = inverse_weight
      )
    )
    names(directed_distance) <- igraph::V(cluster_graph)$name

    undirected_graph <- igraph::as_undirected(cluster_graph, mode = "collapse")
    undirected_inverse_weight <- 1 / pmax(igraph::E(undirected_graph)$weight, eps)
    undirected_distance <- as.numeric(
      igraph::distances(
        undirected_graph,
        v = root_cluster,
        to = igraph::V(undirected_graph),
        mode = "all",
        weights = undirected_inverse_weight
      )
    )
    names(undirected_distance) <- igraph::V(undirected_graph)$name
  }

  cluster_distance <- directed_distance
  use_fallback <- !is.finite(cluster_distance) & is.finite(undirected_distance)
  cluster_distance[use_fallback] <- undirected_distance[use_fallback]
  cluster_pseudotime <- ifelse(is.finite(cluster_distance), cluster_distance, NA_real_)
  finite_pt <- sort(unique(cluster_pseudotime[!is.na(cluster_pseudotime)]))
  positive_diffs <- diff(finite_pt)
  positive_diffs <- positive_diffs[positive_diffs > 0]
  local_step <- if (length(positive_diffs) > 0L) min(positive_diffs) else 1

  cell_pseudotime <- rep(NA_real_, length(clusters))
  names(cell_pseudotime) <- rownames(meta)

  for (cluster_name in levels(clusters)) {
    idx <- which(as.character(clusters) == cluster_name)
    base_pt <- cluster_pseudotime[cluster_name]
    if (length(idx) == 0L || is.na(base_pt)) {
      next
    }

    time_in_cluster <- time_numeric[idx]
    if (length(unique(time_in_cluster)) <= 1L) {
      within_rank <- rep(0, length(idx))
    } else {
      within_rank <- (rank(time_in_cluster, ties.method = "average") - 1) / (length(idx) - 1)
    }

    cell_pseudotime[idx] <- base_pt + within_cluster_fraction * local_step * within_rank
  }

  cluster_pseudotime_scaled <- .scale_to_unit(cluster_pseudotime)
  names(cluster_pseudotime_scaled) <- names(cluster_pseudotime)
  cell_pseudotime_scaled <- .scale_to_unit(cell_pseudotime)
  names(cell_pseudotime_scaled) <- names(cell_pseudotime)

  cluster_pt_table <- data.frame(
    cluster = names(cluster_pseudotime),
    pseudotime = as.numeric(cluster_pseudotime),
    scaled_pseudotime = as.numeric(cluster_pseudotime_scaled[names(cluster_pseudotime)]),
    median_time = as.numeric(cluster_median_time[names(cluster_pseudotime)]),
    size = as.integer(cluster_size[names(cluster_pseudotime)]),
    reachable = !is.na(cluster_pseudotime),
    reachable_directed = is.finite(directed_distance[names(cluster_pseudotime)]),
    used_undirected_fallback = as.logical(use_fallback[names(cluster_pseudotime)]),
    is_root = names(cluster_pseudotime) == root_cluster,
    stringsAsFactors = FALSE
  )

  cluster_pt_table <- cluster_pt_table[order(cluster_pt_table$pseudotime, na.last = TRUE), , drop = FALSE]

  list(
    root_cluster = root_cluster,
    cluster_pseudotime = cluster_pt_table,
    cell_pseudotime = cell_pseudotime,
    cell_pseudotime_scaled = cell_pseudotime_scaled,
    within_cluster_step = local_step
  )
}
