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
