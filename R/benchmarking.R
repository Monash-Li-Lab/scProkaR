#' Benchmark integration quality with scIB-inspired metrics
#'
#' Computes a compact R-native metric panel across available reduced dimensions
#' and summarizes the results with scores and optional plots.
#'
#' @param sce A `SingleCellExperiment`.
#' @param batch_col Column in `colData(sce)` containing batch labels.
#' @param label_col Optional biological label column for conservation metrics.
#' @param methods Optional vector of integration method names or reducedDim names
#'   to benchmark. When `NULL`, all stored integration results are used.
#' @param metrics Metric names to calculate.
#' @param return_plots If `TRUE`, return summary ggplot objects.
#'
#' Any embedding already present in `reducedDims(sce)` can be benchmarked by
#' passing its reduced-dimension name in `methods`, and external embeddings can
#' be registered with [RegisterIntegrationEmbedding()].
#'
#' @return A list with `scores`, `ranking`, and optional `plots`.
#' @export
#' @examples
#' set.seed(1)
#' sce <- simulate_tata_multidrug_sce(
#'     n_cells = 300, n_features = 60, n_pcs = 10)
#' pca <- SingleCellExperiment::reducedDim(sce, "PCA")
#'
#' # Use the drug condition as batch and a quick clustering as the
#' # biological label to conserve.
#' sce$batch <- sce$condition
#' sce$cell_type <- factor(stats::kmeans(pca, centers = 3)$cluster)
#'
#' # Register a second embedding so it is scored alongside the PCA.
#' sce <- RegisterIntegrationEmbedding(
#'     sce, pca[, seq_len(5)], method_name = "pca_top5")
#'
#' bench <- BenchmarkIntegration(
#'     sce,
#'     batch_col = "batch",
#'     label_col = "cell_type",
#'     methods = c("PCA", "pca_top5"),
#'     return_plots = FALSE)
#' head(bench$scores)
#' bench$ranking
BenchmarkIntegration <- function(
    sce,
    batch_col,
    label_col = NULL,
    methods = NULL,
    metrics = c(
      "pcr_batch", "ilisi", "clisi", "kbet_like",
      "graph_connectivity", "ari", "nmi", "silhouette"
    ),
    return_plots = TRUE
) {
  .scprokar_match_columns(sce, batch_col, label = "batch")
  if (!is.null(label_col)) {
    .scprokar_match_columns(sce, label_col, label = "label")
  }

  reductions <- .scprokar_reduction_lookup(sce, methods = methods)
  if (length(reductions) == 0) {
    stop("No reduced dimensions are available to benchmark.", call. = FALSE)
  }

  batch <- as.factor(SummarizedExperiment::colData(sce)[[batch_col]])
  labels <- if (!is.null(label_col)) as.factor(SummarizedExperiment::colData(sce)[[label_col]]) else NULL
  scores <- do.call(
    rbind,
    lapply(names(reductions), function(method_name) {
      reduction_name <- reductions[[method_name]]
      emb <- SingleCellExperiment::reducedDim(sce, reduction_name)
      .scprokar_score_reduction(
        emb = emb,
        batch = batch,
        labels = labels,
        method = method_name,
        reduction = reduction_name,
        metrics = metrics
      )
    })
  )

  ranking <- .scprokar_rank_scores(scores)
  plots <- NULL
  if (isTRUE(return_plots)) {
    plots <- .scprokar_benchmark_plots(scores, ranking)
  }

  meta <- .scprokar_get_metadata(sce)
  meta$benchmark <- list(
    scores = scores,
    ranking = ranking,
    batch_col = batch_col,
    label_col = label_col
  )
  sce <- .scprokar_set_metadata(sce, meta)

  list(sce = sce, scores = scores, ranking = ranking, plots = plots)
}

#' @keywords internal
.scprokar_score_reduction <- function(emb, batch, labels, method, reduction, metrics) {
  emb <- as.matrix(emb)
  k <- min(15, nrow(emb) - 1)
  idx <- .scprokar_neighbor_index(emb, k = k)
  graph <- .scprokar_build_knn_graph(emb, k = k)
  cluster_assign <- .scprokar_cluster_for_metrics(emb, labels)

  metric_values <- lapply(metrics, function(metric) {
    out <- switch(
      metric,
      pcr_batch = .scprokar_metric_pcr_batch(emb, batch),
      ilisi = .scprokar_metric_lisi(idx, batch, mode = "mixing"),
      clisi = .scprokar_metric_lisi(idx, labels, mode = "purity"),
      kbet_like = .scprokar_metric_kbet_like(idx, batch),
      graph_connectivity = .scprokar_metric_graph_connectivity(graph, labels),
      ari = .scprokar_metric_ari(cluster_assign, labels),
      nmi = .scprokar_metric_nmi(cluster_assign, labels),
      silhouette = .scprokar_metric_silhouette(emb, labels),
      stop("Unknown metric: ", metric, call. = FALSE)
    )
    data.frame(
      method = method,
      reduction = reduction,
      metric = metric,
      value = out$value,
      score = out$score,
      higher_better = out$higher_better,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, metric_values)
}

#' @keywords internal
.scprokar_cluster_for_metrics <- function(emb, labels) {
  if (is.null(labels) || length(unique(labels)) < 2) {
    return(factor(rep("cluster_1", nrow(emb))))
  }
  k <- length(unique(labels))
  factor(stats::kmeans(emb, centers = k, iter.max = 25)$cluster)
}

#' @keywords internal
.scprokar_metric_pcr_batch <- function(emb, batch) {
  if (length(unique(batch)) < 2) {
    return(list(value = 1, score = 1, higher_better = TRUE))
  }
  per_dim <- apply(emb, 2, function(dim_values) {
    fit <- stats::lm(dim_values ~ batch)
    summary(fit)$r.squared
  })
  score <- 1 - mean(per_dim)
  list(value = score, score = score, higher_better = TRUE)
}

#' @keywords internal
.scprokar_metric_lisi <- function(idx, labels, mode = c("mixing", "purity")) {
  mode <- match.arg(mode)
  if (is.null(labels) || ncol(idx) == 0) {
    return(list(value = NA_real_, score = NA_real_, higher_better = TRUE))
  }
  labels <- as.factor(labels)
  n_levels <- length(levels(labels))
  if (n_levels <= 1) {
    return(list(value = 1, score = 1, higher_better = TRUE))
  }
  lisi <- apply(idx, 1, function(neighbors) {
    props <- prop.table(table(labels[neighbors]))
    1 / sum(props^2)
  })
  if (mode == "mixing") {
    score <- mean((lisi - 1) / (n_levels - 1))
    return(list(value = mean(lisi), score = score, higher_better = TRUE))
  }
  score <- mean(1 / lisi)
  list(value = mean(lisi), score = score, higher_better = TRUE)
}

#' @keywords internal
.scprokar_metric_kbet_like <- function(idx, batch) {
  if (ncol(idx) == 0 || length(unique(batch)) < 2) {
    return(list(value = 1, score = 1, higher_better = TRUE))
  }
  batch <- as.factor(batch)
  global <- prop.table(table(batch))
  deviation <- apply(idx, 1, function(neighbors) {
    local <- prop.table(table(factor(batch[neighbors], levels = levels(batch))))
    sum(abs(local - global)) / 2
  })
  score <- 1 - mean(deviation)
  list(value = score, score = score, higher_better = TRUE)
}

#' @keywords internal
.scprokar_metric_graph_connectivity <- function(graph, labels) {
  if (is.null(labels)) {
    return(list(value = NA_real_, score = NA_real_, higher_better = TRUE))
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    return(list(value = NA_real_, score = NA_real_, higher_better = TRUE))
  }
  labels <- as.factor(labels)
  g <- igraph::graph_from_adjacency_matrix(graph, mode = "undirected", diag = FALSE)
  connectivity <- vapply(levels(labels), function(label) {
    nodes <- which(labels == label)
    if (length(nodes) <= 1) {
      return(1)
    }
    subgraph <- igraph::induced_subgraph(g, vids = nodes)
    comps <- igraph::components(subgraph)$csize
    max(comps) / length(nodes)
  }, numeric(1))
  score <- mean(connectivity)
  list(value = score, score = score, higher_better = TRUE)
}

#' @keywords internal
.scprokar_metric_ari <- function(cluster_assign, labels) {
  if (is.null(labels)) {
    return(list(value = NA_real_, score = NA_real_, higher_better = TRUE))
  }
  value <- .scprokar_adjusted_rand_index(cluster_assign, labels)
  list(value = value, score = value, higher_better = TRUE)
}

#' @keywords internal
.scprokar_metric_nmi <- function(cluster_assign, labels) {
  if (is.null(labels)) {
    return(list(value = NA_real_, score = NA_real_, higher_better = TRUE))
  }
  value <- .scprokar_nmi(cluster_assign, labels)
  list(value = value, score = value, higher_better = TRUE)
}

#' @keywords internal
.scprokar_metric_silhouette <- function(emb, labels) {
  if (is.null(labels) || length(unique(labels)) < 2) {
    return(list(value = NA_real_, score = NA_real_, higher_better = TRUE))
  }
  if (!requireNamespace("cluster", quietly = TRUE)) {
    return(list(value = NA_real_, score = NA_real_, higher_better = TRUE))
  }
  dist_obj <- stats::dist(emb)
  sil <- cluster::silhouette(as.integer(as.factor(labels)), dist_obj)
  avg <- mean(sil[, "sil_width"])
  list(value = avg, score = (avg + 1) / 2, higher_better = TRUE)
}

#' @keywords internal
.scprokar_rank_scores <- function(scores) {
  method_scores <- stats::aggregate(score ~ method, data = scores, FUN = function(x) mean(x, na.rm = TRUE))
  method_scores <- method_scores[order(method_scores$score, decreasing = TRUE), , drop = FALSE]
  method_scores$rank <- seq_len(nrow(method_scores))
  rownames(method_scores) <- NULL
  method_scores
}

#' @keywords internal
.scprokar_benchmark_plots <- function(scores, ranking) {
  scores$category <- .scprokar_metric_category(scores$metric)
  category_scores <- stats::aggregate(
    score ~ method + category,
    data = scores,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(category_scores)[names(category_scores) == "score"] <- "category_score"

  batch_scores <- category_scores[category_scores$category == "batch_removal", c("method", "category_score"), drop = FALSE]
  bio_scores <- category_scores[category_scores$category == "bio_conservation", c("method", "category_score"), drop = FALSE]
  names(batch_scores)[2] <- "batch_removal"
  names(bio_scores)[2] <- "bio_conservation"
  tradeoff <- merge(batch_scores, bio_scores, by = "method", all = TRUE)

  list(
    overall = ggplot2::ggplot(ranking, ggplot2::aes(x = stats::reorder(method, score), y = score, fill = method)) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Method", y = "Mean score", title = "Integration ranking") +
      ggplot2::theme_minimal(),
    key_metrics = ggplot2::ggplot(scores, ggplot2::aes(x = method, y = score, fill = method)) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::facet_wrap(~metric, scales = "free_y") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Method", y = "Metric score", title = "Key integration metric comparison") +
      ggplot2::theme_minimal(),
    batch_removal = ggplot2::ggplot(
      category_scores[category_scores$category == "batch_removal", , drop = FALSE],
      ggplot2::aes(x = stats::reorder(method, category_score), y = category_score, fill = method)
    ) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Method", y = "Mean batch-removal score", title = "Batch-removal comparison") +
      ggplot2::theme_minimal(),
    bio_conservation = ggplot2::ggplot(
      category_scores[category_scores$category == "bio_conservation", , drop = FALSE],
      ggplot2::aes(x = stats::reorder(method, category_score), y = category_score, fill = method)
    ) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Method", y = "Mean bio-conservation score", title = "Bio-conservation comparison") +
      ggplot2::theme_minimal(),
    tradeoff = ggplot2::ggplot(
      tradeoff,
      ggplot2::aes(x = batch_removal, y = bio_conservation, label = method, color = method)
    ) +
      ggplot2::geom_point(size = 3, show.legend = FALSE) +
      ggplot2::geom_text(vjust = -0.7, show.legend = FALSE) +
      ggplot2::labs(
        x = "Batch-removal score",
        y = "Bio-conservation score",
        title = "Batch-removal versus bio-conservation"
      ) +
      ggplot2::theme_minimal(),
    heatmap = ggplot2::ggplot(scores, ggplot2::aes(x = metric, y = method, fill = score)) +
      ggplot2::geom_tile() +
      ggplot2::labs(x = "Metric", y = "Method", title = "Metric score heatmap") +
      ggplot2::theme_minimal()
  )
}

#' @keywords internal
.scprokar_metric_category <- function(metric) {
  batch_metrics <- c("pcr_batch", "ilisi", "kbet_like", "graph_connectivity")
  ifelse(metric %in% batch_metrics, "batch_removal", "bio_conservation")
}
