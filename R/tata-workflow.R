#' Run the full TATA workflow
#'
#' Run Time Aware Trajectory Analysis (TATA) end-to-end on a
#' `SingleCellExperiment`.
#'
#' @param sce A `SingleCellExperiment`.
#' @param dimred Reduced dimension to use for graph construction.
#' @param dims Optional subset of dimensions.
#' @param time_col Column containing experimental time labels.
#' @param use_time Logical indicating whether known experimental time should be
#'   used during graph abstraction. Set `FALSE` to run a topology-only TATA
#'   variant.
#' @param cluster_col Name of the column where inferred cluster labels will be
#'   stored.
#' @param k Number of neighbours for the cell-cell kNN graph.
#' @param cluster_method Graph clustering method. Choices are:
#' - `"louvain"`: fast modularity-based clustering and the default option.
#' - `"leiden"`: Leiden community detection if available in the installed
#'   `igraph`; otherwise TATA falls back to Louvain.
#' - `"walktrap"`: random-walk clustering that can give slightly coarser states.
#' @param alpha Exponent for topology weights.
#' @param beta Exponent for time weights.
#' @param direction_threshold Threshold for direction assignment.
#' @param prune_threshold Minimum TATA weight to retain an abstract graph edge.
#' @param time_weight_mode How the time signal contributes to the final edge
#'   score. Choices are:
#' - `"directional_confidence"`: use direction strength regardless of whether
#'   the label ordering is `a -> b` or `b -> a`; this is the recommended option
#'   for multi-branch graphs.
#' - `"ordered"`: use the original ordered-pair score `(flow_score + 1) / 2`,
#'   which is stricter but more sensitive to the arbitrary ordering of cluster
#'   labels.
#' @param do_pseudotime Logical indicating whether to compute pseudotime.
#' @param root_cluster Optional root cluster for pseudotime.
#' @param refine_clusters Logical indicating whether to apply a conservative
#'   temporal refinement step after graph clustering.
#' @param refine_min_cells Minimum cluster size considered for temporal
#'   refinement.
#' @param refine_min_split_fraction Minimum fraction of cells required in both
#'   temporal subclusters for a split to be accepted.
#' @param do_branch_probs Logical indicating whether to compute terminal-branch
#'   probabilities and entropy/plasticity scores.
#' @param expected_branches Optional expected number of terminal branches. If
#'   supplied, TATA uses it to guide terminal-state selection for downstream
#'   branch probabilities and branch assignment.
#' @param n_pcs_if_missing Number of PCs to compute if `dimred` is missing.
#'
#' @return A list containing the updated `sce`, cell graph, adjacency matrix,
#'   cluster graph, cluster edge table, cluster pseudotime, cell pseudotime, and
#'   parameters used.
#' @export
#'
#' @examples
#' sce <- simulate_tata_sce(n_cells = 1000, n_features = 60)
#' tata <- run_tata(sce, dimred = "PCA", k = 12)
#' names(tata)
run_tata <- function(
    sce,
    dimred = "PCA",
    dims = NULL,
    time_col = "timepoint",
    use_time = TRUE,
    cluster_col = "cluster",
    k = 15,
    cluster_method = c("louvain", "leiden", "walktrap"),
    alpha = 1,
    beta = 1,
    direction_threshold = 0.1,
    prune_threshold = 0.01,
    time_weight_mode = c("directional_confidence", "ordered"),
    do_pseudotime = TRUE,
    root_cluster = NULL,
    refine_clusters = TRUE,
    refine_min_cells = NULL,
    refine_min_split_fraction = 0.20,
    do_branch_probs = TRUE,
    expected_branches = NULL,
    n_pcs_if_missing = 20
) {
    if (!methods::is(sce, "SingleCellExperiment")) {
        stop("`sce` must be a SingleCellExperiment.", call. = FALSE)
    }

    cluster_method <- match.arg(cluster_method)
    time_weight_mode <- match.arg(time_weight_mode)

    if (isTRUE(use_time) && !time_col %in% colnames(SummarizedExperiment::colData(sce))) {
        stop("Time column `", time_col, "` is not present in colData(sce).", call. = FALSE)
    }

    if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0) {
        stop("`alpha` must be a single non-negative numeric value.", call. = FALSE)
    }

    if (!is.numeric(beta) || length(beta) != 1L || beta < 0) {
        stop("`beta` must be a single non-negative numeric value.", call. = FALSE)
    }

    sce <- .ensure_reduced_dim(sce, dimred = dimred, n_pcs = n_pcs_if_missing)
    knn_result <- build_knn_graph(x = sce, dimred = dimred, dims = dims, k = k)

    inferred_cluster <- cluster_graph_states(
        cell_graph = knn_result$graph,
        method = cluster_method
    )

    if (isTRUE(refine_clusters) && isTRUE(use_time)) {
        if (is.null(refine_min_cells)) {
            refine_min_cells <- max(40L, 2L * knn_result$k)
        }

        inferred_cluster <- .refine_clusters_temporally(
            clusters = inferred_cluster,
            timepoint = SummarizedExperiment::colData(sce)[[time_col]],
            embedding = knn_result$embedding,
            min_cluster_size = refine_min_cells,
            min_split_fraction = refine_min_split_fraction
        )
    }

    SummarizedExperiment::colData(sce)[[cluster_col]] <- inferred_cluster

    topology_table <- compute_topology_weights(
        adjacency = knn_result$adjacency,
        clusters = inferred_cluster
    )

    time_table <- if (isTRUE(use_time)) {
        compute_time_weights(
            adjacency = knn_result$adjacency,
            clusters = inferred_cluster,
            timepoint = SummarizedExperiment::colData(sce)[[time_col]]
        )
    } else {
        .make_neutral_time_table(inferred_cluster)
    }

    tata_graph_result <- build_tata_graph(
        topology_table = topology_table,
        time_table = time_table,
        clusters = inferred_cluster,
        embedding = knn_result$embedding,
        timepoint = if (isTRUE(use_time)) SummarizedExperiment::colData(sce)[[time_col]] else rep(NA_real_, ncol(sce)),
        alpha = alpha,
        beta = beta,
        direction_threshold = direction_threshold,
        prune_threshold = prune_threshold,
        time_weight_mode = time_weight_mode
    )

    if (isTRUE(do_pseudotime)) {
        pseudotime_result <- compute_tata_pseudotime(
            cluster_graph = tata_graph_result$graph,
            sce = sce,
            cluster_col = cluster_col,
            time_col = if (isTRUE(use_time)) time_col else NULL,
            root_cluster = root_cluster
        )

        SummarizedExperiment::colData(sce)$tata_pseudotime <- pseudotime_result$cell_pseudotime
        SummarizedExperiment::colData(sce)$tata_pseudotime_scaled <- pseudotime_result$cell_pseudotime_scaled
    } else {
        pseudotime_result <- list(
            root_cluster = NULL,
            cluster_pseudotime = NULL,
            cell_pseudotime = NULL,
            cell_pseudotime_scaled = NULL,
            within_cluster_step = NULL
        )
    }

    if (isTRUE(do_branch_probs)) {
        branch_result <- compute_tata_branch_probabilities(
            tata_result = list(
                sce = sce,
                cluster_graph = tata_graph_result$graph,
                cluster_pseudotime = pseudotime_result$cluster_pseudotime,
                parameters = list(
                    cluster_col = cluster_col,
                    root_cluster = pseudotime_result$root_cluster
                )
            ),
            expected_branches = expected_branches
        )

        prob_cols <- branch_result$branch_probability_columns
        prob_values <- branch_result$branch_probabilities[, prob_cols, drop = FALSE]
        rownames(prob_values) <- branch_result$branch_probabilities$cell_id
        prob_values <- prob_values[colnames(sce), , drop = FALSE]

        for (col_name in prob_cols) {
            SummarizedExperiment::colData(sce)[[col_name]] <- prob_values[[col_name]]
        }
        SummarizedExperiment::colData(sce)$tata_branch_entropy <- branch_result$branch_probabilities$tata_branch_entropy
        SummarizedExperiment::colData(sce)$tata_branch_plasticity <- branch_result$branch_probabilities$tata_branch_plasticity
        SummarizedExperiment::colData(sce)$tata_branch_commitment <- branch_result$branch_probabilities$tata_branch_commitment
        SummarizedExperiment::colData(sce)$tata_max_branch_probability <- branch_result$branch_probabilities$tata_max_branch_probability
        SummarizedExperiment::colData(sce)$tata_branch_assignment <- branch_result$branch_probabilities$tata_branch_assignment
    } else {
        branch_result <- list(
            terminal_clusters = NULL,
            transition_matrix = NULL,
            cluster_branch_probabilities = NULL,
            branch_probabilities = NULL,
            branch_probability_columns = character(0)
        )
    }

    cell_space <- .compute_tata_cell_space(
        embedding = knn_result$embedding,
        pseudotime = pseudotime_result$cell_pseudotime_scaled,
        branch_probabilities = branch_result$branch_probabilities,
        branch_probability_columns = branch_result$branch_probability_columns
    )
    SingleCellExperiment::reducedDim(sce, "TATA") <- cell_space

    parameters <- list(
        dimred = dimred,
        dims = dims,
        time_col = time_col,
        use_time = use_time,
        cluster_col = cluster_col,
        k = knn_result$k,
        knn_engine = knn_result$knn_engine,
        cluster_method = cluster_method,
        alpha = alpha,
        beta = beta,
        direction_threshold = direction_threshold,
        prune_threshold = prune_threshold,
        time_weight_mode = time_weight_mode,
        do_pseudotime = do_pseudotime,
        refine_clusters = refine_clusters,
        refine_min_cells = refine_min_cells,
        refine_min_split_fraction = refine_min_split_fraction,
        do_branch_probs = do_branch_probs,
        expected_branches = expected_branches,
        root_cluster = pseudotime_result$root_cluster
    )

    S4Vectors::metadata(sce)$TATA <- list(
        parameters = parameters,
        cluster_edge_table = tata_graph_result$edge_table,
        cluster_vertex_table = tata_graph_result$vertex_table,
        cluster_pseudotime = pseudotime_result$cluster_pseudotime,
        terminal_clusters = branch_result$terminal_clusters,
        cluster_branch_probabilities = branch_result$cluster_branch_probabilities
    )

    list(
        sce = sce,
        cell_graph = knn_result$graph,
        adjacency_matrix = knn_result$adjacency,
        cluster_graph = tata_graph_result$graph,
        cluster_edge_table = tata_graph_result$edge_table,
        cluster_vertex_table = tata_graph_result$vertex_table,
        cluster_pseudotime = pseudotime_result$cluster_pseudotime,
        cell_pseudotime = pseudotime_result$cell_pseudotime,
        cell_pseudotime_scaled = pseudotime_result$cell_pseudotime_scaled,
        cell_space = cell_space,
        terminal_clusters = branch_result$terminal_clusters,
        cluster_branch_probabilities = branch_result$cluster_branch_probabilities,
        branch_probabilities = branch_result$branch_probabilities,
        transition_matrix = branch_result$transition_matrix,
        parameters = parameters,
        topology_table = topology_table,
        time_table = time_table
    )
}
