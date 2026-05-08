# Internal helpers used across the TATA workflow.

.coerce_time_to_numeric <- function(timepoint) {
  if (inherits(timepoint, c("Date", "POSIXct", "POSIXt"))) {
    return(as.numeric(timepoint))
  }

  if (is.numeric(timepoint) || is.integer(timepoint)) {
    return(as.numeric(timepoint))
  }

  if (is.factor(timepoint)) {
    suppressWarnings(numeric_levels <- as.numeric(as.character(timepoint)))
    if (!anyNA(numeric_levels)) {
      return(numeric_levels)
    }
    return(as.numeric(timepoint))
  }

  if (is.character(timepoint)) {
    suppressWarnings(numeric_values <- as.numeric(timepoint))
    if (!anyNA(numeric_values)) {
      return(numeric_values)
    }
    ordered_levels <- sort(unique(timepoint))
    return(match(timepoint, ordered_levels))
  }

  stop("`timepoint` could not be converted to numeric.", call. = FALSE)
}


.ensure_reduced_dim <- function(sce, dimred = "PCA", n_pcs = 20) {
  if (!methods::is(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment.", call. = FALSE)
  }

  if (dimred %in% SingleCellExperiment::reducedDimNames(sce)) {
    return(sce)
  }

  assay_name <- NULL
  if ("logcounts" %in% SummarizedExperiment::assayNames(sce)) {
    assay_name <- "logcounts"
  } else if ("counts" %in% SummarizedExperiment::assayNames(sce)) {
    assay_name <- "counts"
  }

  if (is.null(assay_name)) {
    stop(
      "No reduced dimension named `", dimred, "` was found, and neither ",
      "`logcounts` nor `counts` assays are available to compute PCA.",
      call. = FALSE
    )
  }

  mat <- t(as.matrix(SummarizedExperiment::assay(sce, assay_name)))
  pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
  n_keep <- min(n_pcs, ncol(pca$x))
  pcs <- pca$x[, seq_len(n_keep), drop = FALSE]
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  SingleCellExperiment::reducedDim(sce, dimred) <- pcs
  sce
}


.get_embedding_matrix <- function(x, dimred = "PCA", dims = NULL) {
  if (methods::is(x, "SingleCellExperiment")) {
    if (!dimred %in% SingleCellExperiment::reducedDimNames(x)) {
      stop("Reduced dimension `", dimred, "` is not present in `sce`.", call. = FALSE)
    }
    embedding <- as.matrix(SingleCellExperiment::reducedDim(x, dimred))
    rownames(embedding) <- colnames(x)
  } else {
    embedding <- as.matrix(x)
  }

  if (!is.null(dims)) {
    embedding <- embedding[, dims, drop = FALSE]
  }

  if (nrow(embedding) < 2L) {
    stop("The embedding must contain at least two cells.", call. = FALSE)
  }

  if (is.null(rownames(embedding))) {
    rownames(embedding) <- paste0("cell_", seq_len(nrow(embedding)))
  }

  embedding
}


.compute_knn <- function(embedding, k = 15) {
  n_cells <- nrow(embedding)
  k <- max(1, min(as.integer(k), n_cells - 1L))

  if (requireNamespace("FNN", quietly = TRUE)) {
    knn <- FNN::get.knn(data = embedding, k = k)
    return(list(
      nn.index = knn$nn.index,
      nn.dist = knn$nn.dist,
      engine = "FNN"
    ))
  }

  if (requireNamespace("BiocNeighbors", quietly = TRUE)) {
    raw_knn <- BiocNeighbors::findKNN(embedding, k = k + 1L)
    nn.index <- matrix(NA_integer_, nrow = n_cells, ncol = k)
    nn.dist <- matrix(NA_real_, nrow = n_cells, ncol = k)

    for (i in seq_len(n_cells)) {
      keep <- raw_knn$index[i, ] != i
      idx_i <- raw_knn$index[i, keep]
      dist_i <- raw_knn$distance[i, keep]
      if (length(idx_i) < k) {
        stop("BiocNeighbors did not return enough non-self neighbors.", call. = FALSE)
      }
      nn.index[i, ] <- idx_i[seq_len(k)]
      nn.dist[i, ] <- dist_i[seq_len(k)]
    }

    return(list(
      nn.index = nn.index,
      nn.dist = nn.dist,
      engine = "BiocNeighbors"
    ))
  }

  distance_matrix <- as.matrix(stats::dist(embedding))
  diag(distance_matrix) <- Inf
  nn.index <- t(apply(distance_matrix, 1L, order))[ , seq_len(k), drop = FALSE]
  nn.dist <- matrix(
    distance_matrix[cbind(rep(seq_len(n_cells), each = k), as.vector(t(nn.index)))],
    nrow = n_cells,
    ncol = k,
    byrow = TRUE
  )

  list(
    nn.index = nn.index,
    nn.dist = nn.dist,
    engine = "base::dist"
  )
}


.compute_cluster_centroids <- function(embedding, clusters) {
  xy <- embedding[, seq_len(min(2L, ncol(embedding))), drop = FALSE]
  if (ncol(xy) == 1L) {
    xy <- cbind(xy, 0)
  }
  colnames(xy) <- c("x", "y")

  stats::aggregate(
    xy,
    by = list(cluster = as.character(clusters)),
    FUN = stats::median
  )
}


.scale_to_unit <- function(x) {
  finite_idx <- which(!is.na(x))
  scaled <- rep(NA_real_, length(x))

  if (length(finite_idx) == 0L) {
    return(scaled)
  }

  if (length(unique(x[finite_idx])) == 1L) {
    scaled[finite_idx] <- 0
    return(scaled)
  }

  x_range <- range(x[finite_idx])
  scaled[finite_idx] <- (x[finite_idx] - x_range[1]) / diff(x_range)
  scaled
}


.typical_time_step <- function(timepoint, fallback = 1) {
  time_numeric <- sort(unique(.coerce_time_to_numeric(timepoint)))
  time_numeric <- time_numeric[is.finite(time_numeric)]
  diffs <- diff(time_numeric)
  diffs <- diffs[diffs > 0]

  if (length(diffs) == 0L) {
    return(as.numeric(fallback))
  }

  stats::median(diffs)
}


.time_distribution_overlap <- function(time_a, time_b) {
  time_a <- .coerce_time_to_numeric(time_a)
  time_b <- .coerce_time_to_numeric(time_b)
  time_a <- time_a[is.finite(time_a)]
  time_b <- time_b[is.finite(time_b)]

  if (length(time_a) == 0L || length(time_b) == 0L) {
    return(1)
  }

  support <- sort(unique(c(time_a, time_b)))
  prop_a <- table(factor(time_a, levels = support)) / length(time_a)
  prop_b <- table(factor(time_b, levels = support)) / length(time_b)
  sum(pmin(as.numeric(prop_a), as.numeric(prop_b)))
}


.refine_clusters_temporally <- function(
    clusters,
    timepoint,
    embedding = NULL,
    min_cluster_size = 60L,
    min_split_fraction = 0.20,
    min_center_gap = 1.1,
    min_between_ratio = 0.30,
    min_embedding_gap = 0.55,
    seed = 1L) {
  time_numeric <- .coerce_time_to_numeric(timepoint)
  if (!all(is.finite(time_numeric))) {
    return(factor(clusters))
  }

  typical_step <- .typical_time_step(time_numeric)
  clusters <- as.character(clusters)
  refined <- clusters

  for (cluster_name in unique(clusters)) {
    idx <- which(clusters == cluster_name)
    n_cluster <- length(idx)

    if (n_cluster < max(6L, as.integer(min_cluster_size))) {
      next
    }

    cluster_time <- time_numeric[idx]
    if (stats::IQR(cluster_time) < typical_step) {
      next
    }

    set.seed(seed + match(cluster_name, unique(clusters)))
    split_fit <- tryCatch(
      stats::kmeans(cluster_time, centers = 2L, nstart = 20L),
      error = function(e) NULL
    )

    if (is.null(split_fit)) {
      next
    }

    split_sizes <- tabulate(split_fit$cluster, nbins = 2L)
    if (min(split_sizes) < min_split_fraction * n_cluster) {
      next
    }

    ordered_centers <- order(as.numeric(split_fit$centers))
    center_gap <- diff(sort(as.numeric(split_fit$centers)))
    if (length(center_gap) == 0L || center_gap < min_center_gap * typical_step) {
      next
    }

    between_ratio <- split_fit$betweenss / max(split_fit$totss, 1e-8)
    if (!is.finite(between_ratio) || between_ratio < min_between_ratio) {
      next
    }

    if (!is.null(embedding)) {
      emb_cluster <- embedding[idx, , drop = FALSE]
      emb_cluster <- emb_cluster[, seq_len(min(5L, ncol(emb_cluster))), drop = FALSE]
      if (ncol(emb_cluster) > 0L) {
        emb_cluster <- scale(emb_cluster)
        split_groups <- match(split_fit$cluster, ordered_centers)
        centroid_1 <- colMeans(emb_cluster[split_groups == 1L, , drop = FALSE], na.rm = TRUE)
        centroid_2 <- colMeans(emb_cluster[split_groups == 2L, , drop = FALSE], na.rm = TRUE)
        embedding_gap <- sqrt(sum((centroid_1 - centroid_2) ^ 2))
        if (!is.finite(embedding_gap) || embedding_gap < min_embedding_gap) {
          next
        }
      }
    }

    ordered_labels <- match(split_fit$cluster, ordered_centers)
    refined[idx] <- paste0(cluster_name, "_", ordered_labels)
  }

  refined_levels <- unique(refined)
  median_time <- tapply(time_numeric, refined, stats::median)
  refined_levels <- names(sort(median_time[refined_levels], na.last = TRUE))
  relabel_map <- stats::setNames(
    paste0("C", seq_along(refined_levels)),
    refined_levels
  )

  factor(relabel_map[refined], levels = unname(relabel_map))
}


.compute_tata_cell_space <- function(
    embedding,
    pseudotime = NULL,
    branch_probabilities = NULL,
    branch_probability_columns = character(0),
    n_components = 5L) {
  .safe_scale <- function(mat) {
    mat <- as.matrix(mat)
    if (!is.numeric(mat)) {
      storage.mode(mat) <- "numeric"
    }
    keep <- apply(mat, 2L, function(x) {
      finite_x <- x[is.finite(x)]
      length(finite_x) > 1L && stats::sd(finite_x) > 0
    })
    if (!any(keep)) {
      return(matrix(numeric(0), nrow = nrow(mat), ncol = 0L))
    }
    mat <- mat[, keep, drop = FALSE]
    for (j in seq_len(ncol(mat))) {
      col_j <- mat[, j]
      finite_j <- col_j[is.finite(col_j)]
      fill_value <- if (length(finite_j) > 0L) stats::median(finite_j) else 0
      col_j[!is.finite(col_j)] <- fill_value
      sd_j <- stats::sd(col_j)
      if (!is.finite(sd_j) || sd_j <= 0) {
        col_j[] <- 0
      } else {
        col_j <- (col_j - mean(col_j)) / sd_j
      }
      mat[, j] <- col_j
    }
    mat
  }

  embedding <- as.matrix(embedding)
  if (!is.numeric(embedding)) {
    storage.mode(embedding) <- "numeric"
  }

  keep_dims <- seq_len(min(5L, ncol(embedding)))
  feature_blocks <- list(.safe_scale(embedding[, keep_dims, drop = FALSE]))

  if (!is.null(pseudotime)) {
    feature_blocks[[length(feature_blocks) + 1L]] <- matrix(
      .scale_to_unit(as.numeric(pseudotime)),
      ncol = 1L
    )
  }

  if (!is.null(branch_probabilities) && length(branch_probability_columns) > 0L) {
    prob_mat <- as.matrix(branch_probabilities[, branch_probability_columns, drop = FALSE])
    if (!is.numeric(prob_mat)) {
      storage.mode(prob_mat) <- "numeric"
    }
    prob_mat <- .safe_scale(prob_mat)
    feature_blocks[[length(feature_blocks) + 1L]] <- 0.60 * prob_mat
  }

  feature_mat <- do.call(cbind, feature_blocks)
  if (ncol(feature_mat) == 0L) {
    return(matrix(0, nrow = nrow(embedding), ncol = 1L))
  }
  feature_mat <- .safe_scale(feature_mat)

  if (ncol(feature_mat) == 0L) {
    return(matrix(0, nrow = nrow(embedding), ncol = 1L))
  }

  if (ncol(feature_mat) == 1L) {
    out <- matrix(feature_mat[, 1], ncol = 1L)
    colnames(out) <- "TATA1"
    rownames(out) <- rownames(embedding)
    return(out)
  }

  fit <- stats::prcomp(feature_mat, center = TRUE, scale. = FALSE)
  n_keep <- min(as.integer(n_components), ncol(fit$x))
  out <- fit$x[, seq_len(n_keep), drop = FALSE]
  colnames(out) <- paste0("TATA", seq_len(ncol(out)))
  rownames(out) <- rownames(embedding)
  out
}
