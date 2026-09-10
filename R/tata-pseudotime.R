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
#'
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' tata <- run_tata(sce, dimred = "PCA", k = 15)
#'
#' # Recompute pseudotime directly from the abstracted cluster graph.
#' pt <- compute_tata_pseudotime(
#'     cluster_graph = tata$cluster_graph,
#'     sce = tata$sce,
#'     cluster_col = "cluster",
#'     time_col = "timepoint"
#' )
#'
#' pt$root_cluster
#' head(pt$cluster_pseudotime)
#' summary(pt$cell_pseudotime_scaled)
compute_tata_pseudotime <- function(
    cluster_graph,
    sce,
    cluster_col = "cluster",
    time_col = "timepoint",
    root_cluster = NULL,
    within_cluster_fraction = 0.5,
    eps = 1e-8
) {
    if (!inherits(cluster_graph, "igraph")) {
        stop("`cluster_graph` must be an igraph object.", call. = FALSE)
    }

    meta <- as.data.frame(SummarizedExperiment::colData(sce))
    if (!cluster_col %in% colnames(meta)) {
        stop(
            "Cluster column `", cluster_col, "` was not found in colData.",
            call. = FALSE
        )
    }
    if (!is.null(time_col) && !time_col %in% colnames(meta)) {
        stop(
            "Time column `", time_col, "` was not found in colData.",
            call. = FALSE
        )
    }

    clusters <- factor(
        meta[[cluster_col]],
        levels = igraph::V(cluster_graph)$name
    )
    has_time <- !is.null(time_col)
    time_numeric <- if (isTRUE(has_time)) .coerce_time_to_numeric(
        meta[[time_col]]
    ) else rep(NA_real_, nrow(meta))

    cluster_size <- tapply(seq_along(clusters), clusters, length)
    cluster_median_time <- if (isTRUE(has_time)) {
        tapply(time_numeric, clusters, stats::median)
    } else {
        stats::setNames(
            rep(NA_real_, length(levels(clusters))),
            levels(clusters)
        )
    }

    if (is.null(root_cluster)) {
        undirected_graph <- igraph::as_undirected(
            cluster_graph,
            mode = "collapse"
        )
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

        if (isTRUE(has_time)) {
            earliest_time <- min(time_numeric, na.rm = TRUE)
            earliest_fraction <- tapply(
                time_numeric == earliest_time,
                clusters,
                mean
            )
            earliest_fraction[is.na(earliest_fraction)] <- 0

            candidate_mask <- earliest_fraction > 0
            candidates <- names(earliest_fraction)[candidate_mask]
            if (length(candidates) == 0L) {
                earliest_cluster_time <- min(cluster_median_time, na.rm = TRUE)
                candidates <- names(cluster_median_time)[
                    cluster_median_time == earliest_cluster_time
                ]
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
                ), ,
                drop = FALSE
            ]

            root_cluster <- candidate_table$cluster[1]
        } else {
            candidate_table <- data.frame(
                cluster = names(closeness_score),
                closeness = as.numeric(closeness_score),
                size = as.numeric(cluster_size[names(closeness_score)]),
                stringsAsFactors = FALSE
            )
            candidate_table <- candidate_table[
                order(
                    -candidate_table$closeness,
                    -candidate_table$size
                ), ,
                drop = FALSE
            ]
            root_cluster <- candidate_table$cluster[1]
        }
    }

    if (!root_cluster %in% igraph::V(cluster_graph)$name) {
        stop(
            "`root_cluster` is not present in the cluster graph.",
            call. = FALSE
        )
    }

    if (igraph::ecount(cluster_graph) == 0L) {
        directed_distance <- rep(Inf, igraph::vcount(cluster_graph))
        names(directed_distance) <- igraph::V(cluster_graph)$name
        directed_distance[root_cluster] <- 0
        undirected_distance <- directed_distance
        cluster_pseudotime <- directed_distance
    } else if (!isTRUE(has_time)) {
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

        undirected_graph <- igraph::as_undirected(
            cluster_graph,
            mode = "collapse"
        )
        undirected_inverse_weight <- 1 / pmax(
            igraph::E(undirected_graph)$weight,
            eps
        )
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

        cluster_pseudotime <- directed_distance
        use_fallback <- !is.finite(cluster_pseudotime) &
            is.finite(undirected_distance)
        cluster_pseudotime[use_fallback] <- undirected_distance[use_fallback]
    } else {
        edge_df <- igraph::as_data_frame(cluster_graph, what = "edges")
        if (!"ambiguous" %in% colnames(edge_df)) {
            edge_df$ambiguous <- FALSE
        }
        if (!"time_consistency" %in% colnames(edge_df)) {
            edge_df$time_consistency <- 1
        }

        cluster_levels <- igraph::V(cluster_graph)$name
        typical_step <- .typical_time_step(time_numeric)
        median_lookup <- stats::setNames(
            cluster_median_time[cluster_levels],
            cluster_levels
        )
        max_weight <- max(edge_df$weight, na.rm = TRUE)
        if (!is.finite(max_weight) || max_weight <= 0) {
            max_weight <- 1
        }

        order_table <- data.frame(
            cluster = cluster_levels,
            median_time = as.numeric(median_lookup[cluster_levels]),
            stringsAsFactors = FALSE
        )
        order_table$median_time[!is.finite(order_table$median_time)] <- max(
            order_table$median_time[is.finite(order_table$median_time)],
            na.rm = TRUE
        )
        order_table <- order_table[
            order(order_table$median_time, order_table$cluster), ,
            drop = FALSE
        ]
        cluster_order <- unique(
            c(root_cluster, setdiff(order_table$cluster, root_cluster))
        )
        rank_lookup <- stats::setNames(seq_along(cluster_order), cluster_order)

        edge_df$delta_time <- median_lookup[edge_df$to] -
            median_lookup[edge_df$from]
        edge_df$from_rank <- rank_lookup[edge_df$from]
        edge_df$to_rank <- rank_lookup[edge_df$to]
        forward_edge_df <- edge_df[
            edge_df$to_rank > edge_df$from_rank, ,
            drop = FALSE
        ]

        cluster_pseudotime <- stats::setNames(
            rep(NA_real_, length(cluster_levels)),
            cluster_levels
        )
        cluster_pseudotime[root_cluster] <- 0

        time_anchor <- pmax(
            0,
            (median_lookup - median_lookup[root_cluster]) /
                pmax(typical_step, eps)
        )

        for (cluster_name in cluster_order[-1]) {
            incoming <- forward_edge_df[
                forward_edge_df$to == cluster_name, ,
                drop = FALSE
            ]
            if (nrow(incoming) == 0L) {
                next
            }

            pred_pt <- cluster_pseudotime[incoming$from]
            valid <- is.finite(pred_pt)
            if (!any(valid)) {
                next
            }

            incoming <- incoming[valid, , drop = FALSE]
            pred_pt <- pred_pt[valid]

            delta_units <- pmax(
                incoming$delta_time / pmax(typical_step, eps),
                0
            )
            base_increment <- pmax(0.35, delta_units)
            weight_penalty <- 1 + 0.35 * (1 - incoming$weight / max_weight)
            ambiguity_penalty <- ifelse(incoming$ambiguous, 1.35, 1)
            consistency_penalty <- 1 + 0.25 * (
                1 - pmin(pmax(incoming$time_consistency, 0), 1)
            )
            edge_increment <- base_increment * weight_penalty *
                ambiguity_penalty * consistency_penalty

            cluster_pseudotime[cluster_name] <- max(
                pred_pt + edge_increment,
                na.rm = TRUE
            )
        }

        directed_distance <- cluster_pseudotime
        directed_distance[!is.finite(directed_distance)] <- Inf

        undirected_graph <- igraph::as_undirected(
            cluster_graph,
            mode = "collapse"
        )
        undirected_inverse_weight <- 1 / pmax(
            igraph::E(undirected_graph)$weight,
            eps
        )
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

        finite_path <- cluster_pseudotime[is.finite(cluster_pseudotime)]
        max_path <- if (length(finite_path) > 0L) max(finite_path) else 0
        fallback_scale <- .scale_to_unit(undirected_distance)
        fallback_pt <- time_anchor + 0.75 + fallback_scale * (max_path + 0.75)

        use_fallback <- !is.finite(cluster_pseudotime) &
            is.finite(undirected_distance)
        cluster_pseudotime[use_fallback] <- fallback_pt[use_fallback]
        cluster_pseudotime[is.finite(cluster_pseudotime)] <- pmax(
            cluster_pseudotime[is.finite(cluster_pseudotime)],
            0.85 * time_anchor[is.finite(cluster_pseudotime)]
        )
    }

    use_fallback <- !is.finite(directed_distance) &
        is.finite(undirected_distance)
    cluster_pseudotime <- ifelse(
        is.finite(cluster_pseudotime),
        cluster_pseudotime,
        NA_real_
    )
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
        if (
            !isTRUE(has_time) ||
                length(
                    unique(time_in_cluster[is.finite(time_in_cluster)])
                ) <= 1L
        ) {
            within_rank <- rep(0, length(idx))
        } else {
            within_rank <- (
                rank(time_in_cluster, ties.method = "average") - 1
            ) / (length(idx) - 1)
        }

        cell_pseudotime[idx] <- base_pt +
            within_cluster_fraction * local_step * within_rank
    }

    cluster_pseudotime_scaled <- .scale_to_unit(cluster_pseudotime)
    names(cluster_pseudotime_scaled) <- names(cluster_pseudotime)
    cell_pseudotime_scaled <- .scale_to_unit(cell_pseudotime)
    names(cell_pseudotime_scaled) <- names(cell_pseudotime)

    cluster_pt_table <- data.frame(
        cluster = names(cluster_pseudotime),
        pseudotime = as.numeric(cluster_pseudotime),
        scaled_pseudotime = as.numeric(
            cluster_pseudotime_scaled[names(cluster_pseudotime)]
        ),
        median_time = as.numeric(
            cluster_median_time[names(cluster_pseudotime)]
        ),
        size = as.integer(cluster_size[names(cluster_pseudotime)]),
        reachable = !is.na(cluster_pseudotime),
        reachable_directed = is.finite(
            directed_distance[names(cluster_pseudotime)]
        ),
        used_undirected_fallback = as.logical(
            use_fallback[names(cluster_pseudotime)]
        ),
        is_root = names(cluster_pseudotime) == root_cluster,
        stringsAsFactors = FALSE
    )

    cluster_pt_table <- cluster_pt_table[
        order(cluster_pt_table$pseudotime, na.last = TRUE), ,
        drop = FALSE
    ]

    list(
        root_cluster = root_cluster,
        cluster_pseudotime = cluster_pt_table,
        cell_pseudotime = cell_pseudotime,
        cell_pseudotime_scaled = cell_pseudotime_scaled,
        within_cluster_step = local_step
    )
}
