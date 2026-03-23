#' Internal utilities for SCProkaR
#'
#' Helper functions used across object creation, integration, benchmarking, and
#' downstream analysis.
#'
#' @name scprokar-internal
#' @keywords internal
#' @noRd
NULL

#' @keywords internal
.scprokar_require <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message <- paste0("Package '", pkg, "' is required")
    if (!is.null(reason)) {
      message <- paste0(message, " for ", reason)
    }
    stop(paste0(message, ". Install it first."), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.scprokar_package_version <- function() {
  desc <- tryCatch(
    utils::packageDescription("SCProkaR"),
    warning = function(w) NULL,
    error = function(e) NULL
  )
  if (!is.null(desc) && !is.null(desc$Version)) {
    return(as.character(desc$Version))
  }

  desc_file <- system.file("DESCRIPTION", package = "SCProkaR")
  if (!nzchar(desc_file) || !file.exists(desc_file)) {
    desc_file <- file.path(getwd(), "DESCRIPTION")
  }
  if (file.exists(desc_file)) {
    dcf <- tryCatch(read.dcf(desc_file), error = function(e) NULL)
    if (!is.null(dcf) && "Version" %in% colnames(dcf)) {
      return(as.character(dcf[1, "Version"]))
    }
  }

  NA_character_
}

#' @keywords internal
.scprokar_stopifnot_counts <- function(counts) {
  valid <- inherits(counts, "matrix") || inherits(counts, "Matrix")
  if (!valid) {
    stop(
      "`x` must be a matrix-like gene-by-cell count object or a SingleCellExperiment.",
      call. = FALSE
    )
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Counts must have both gene names and cell names.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.scprokar_as_dgC <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  methods::as(x, "dgCMatrix")
}

#' @keywords internal
.scprokar_align_data_frame <- function(x, ids, what) {
  if (is.null(x)) {
    return(data.frame(row.names = ids))
  }
  x <- as.data.frame(x)
  if (nrow(x) != length(ids)) {
    stop("`", what, "` must have ", length(ids), " rows.", call. = FALSE)
  }
  if (is.null(rownames(x))) {
    rownames(x) <- ids
  }
  missing_ids <- setdiff(ids, rownames(x))
  if (length(missing_ids) > 0) {
    stop("`", what, "` is missing identifiers: ", paste(missing_ids, collapse = ", "), call. = FALSE)
  }
  x[ids, , drop = FALSE]
}

#' @keywords internal
.scprokar_get_metadata <- function(sce) {
  meta <- S4Vectors::metadata(sce)$SCProkaR
  if (is.null(meta)) {
    meta <- list()
  }
  meta
}

#' @keywords internal
.scprokar_set_metadata <- function(sce, meta) {
  all_meta <- S4Vectors::metadata(sce)
  all_meta$SCProkaR <- meta
  S4Vectors::metadata(sce) <- all_meta
  sce
}

#' @keywords internal
.scprokar_resolve_reduction_name <- function(sce, reduction) {
  available <- SingleCellExperiment::reducedDimNames(sce)
  if (!length(available)) {
    stop("No reduced dimensions are available in `sce`.", call. = FALSE)
  }
  if (is.null(reduction) || !nzchar(reduction)) {
    stop("`reduction` must be provided.", call. = FALSE)
  }
  if (reduction %in% available) {
    return(reduction)
  }

  lower_available <- tolower(available)
  idx <- which(lower_available == tolower(reduction))
  if (length(idx) == 1L) {
    return(available[[idx]])
  }

  stop(
    "Reduction '", reduction, "' was not found. Available reductions: ",
    paste(available, collapse = ", "),
    call. = FALSE
  )
}

#' @keywords internal
.scprokar_store_step <- function(sce, step, value) {
  meta <- .scprokar_get_metadata(sce)
  meta[[step]] <- value
  .scprokar_set_metadata(sce, meta)
}

#' @keywords internal
.scprokar_match_columns <- function(sce, columns, label) {
  present <- colnames(SummarizedExperiment::colData(sce))
  missing_columns <- setdiff(stats::na.omit(columns), present)
  if (length(missing_columns) > 0) {
    stop(
      "Missing ", label, " column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
.scprokar_fraction_from_features <- function(counts, feature_index) {
  totals <- Matrix::colSums(counts)
  selected <- if (any(feature_index)) {
    Matrix::colSums(counts[feature_index, , drop = FALSE])
  } else {
    rep(0, ncol(counts))
  }
  fraction <- as.numeric(selected / pmax(totals, 1))
  list(
    totals = as.numeric(totals),
    selected = as.numeric(selected),
    fraction = fraction
  )
}

#' @keywords internal
.scprokar_normalize_logcounts <- function(sce, assay_name = "counts", scale_factor = 1e4) {
  counts <- SummarizedExperiment::assay(sce, assay_name)
  libsize <- Matrix::colSums(counts)
  scaling <- scale_factor / pmax(libsize, 1)
  norm <- counts %*% Matrix::Diagonal(x = scaling)
  lognorm <- .scprokar_as_dgC(log1p(norm))
  dimnames(lognorm) <- dimnames(counts)
  SummarizedExperiment::assay(sce, "logcounts", withDimnames = FALSE) <- lognorm
  sce
}

#' @keywords internal
.scprokar_select_features <- function(sce, feature_set = c("hvg", "all"), nfeatures = 2000) {
  feature_set <- match.arg(feature_set)
  if (feature_set == "all") {
    return(rownames(sce))
  }

  current <- SummarizedExperiment::rowData(sce)
  if ("is_hvg" %in% colnames(current) && any(current$is_hvg, na.rm = TRUE)) {
    return(rownames(sce)[which(current$is_hvg)])
  }

  if ("logcounts" %in% SummarizedExperiment::assayNames(sce) &&
      requireNamespace("scran", quietly = TRUE)) {
    var_fit <- scran::modelGeneVar(sce, assay.type = "logcounts")
    ord <- order(var_fit$bio, decreasing = TRUE)
  } else {
    logx <- SummarizedExperiment::assay(sce, "logcounts")
    center <- Matrix::rowMeans(logx)
    variance <- Matrix::rowMeans((logx - center)^2)
    ord <- order(as.numeric(variance), decreasing = TRUE)
  }

  keep <- ord[seq_len(min(length(ord), nfeatures))]
  row_data <- SummarizedExperiment::rowData(sce)
  row_data$is_hvg <- FALSE
  row_data$is_hvg[keep] <- TRUE
  row_data$hvg_rank <- NA_integer_
  row_data$hvg_rank[ord] <- seq_along(ord)
  SummarizedExperiment::rowData(sce) <- row_data
  rownames(sce)[keep]
}

#' @keywords internal
.scprokar_run_pca <- function(
    sce,
    assay_name = "logcounts",
    features = rownames(sce),
    ncomponents = 30,
    reduction_name = "PCA"
) {
  mat <- SummarizedExperiment::assay(sce, assay_name)[features, , drop = FALSE]
  x <- t(as.matrix(mat))
  x <- scale(x, center = TRUE, scale = TRUE)
  x[is.na(x)] <- 0
  rank_k <- min(ncomponents, max(1, ncol(x) - 1))
  full_rank <- min(nrow(x), ncol(x))
  use_irlba <- requireNamespace("irlba", quietly = TRUE) &&
    ncol(x) > 2 &&
    rank_k < full_rank - 1 &&
    rank_k < floor(full_rank * 0.5)

  if (use_irlba) {
    pcs <- irlba::prcomp_irlba(x, n = rank_k, center = FALSE, scale. = FALSE)
  } else {
    pcs <- stats::prcomp(x, rank. = rank_k, center = FALSE, scale. = FALSE)
  }

  emb <- pcs$x
  rownames(emb) <- colnames(sce)
  SingleCellExperiment::reducedDim(sce, reduction_name) <- emb
  sce
}

#' @keywords internal
.scprokar_prepare_reduction <- function(
    sce,
    assay_name = "counts",
    feature_set = c("hvg", "all"),
    dims = 1:30,
    reduction_name = "PCA"
) {
  feature_set <- match.arg(feature_set)
  if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
    sce <- .scprokar_normalize_logcounts(sce, assay_name = assay_name)
  }
  if (!reduction_name %in% SingleCellExperiment::reducedDimNames(sce) ||
      max(dims) > ncol(SingleCellExperiment::reducedDim(sce, reduction_name))) {
    features <- .scprokar_select_features(sce, feature_set = feature_set)
    sce <- .scprokar_run_pca(
      sce,
      assay_name = "logcounts",
      features = features,
      ncomponents = max(dims),
      reduction_name = reduction_name
    )
  }
  sce
}

#' @keywords internal
.scprokar_neighbor_index <- function(x, k = 15) {
  n <- nrow(x)
  if (n <= 1) {
    return(matrix(integer(0), nrow = n, ncol = 0))
  }
  k <- min(k, n - 1)

  if (requireNamespace("RANN", quietly = TRUE)) {
    idx <- RANN::nn2(x, x, k = k + 1)$nn.idx[, -1, drop = FALSE]
    return(idx)
  }

  dmat <- as.matrix(stats::dist(x))
  t(apply(dmat, 1, function(v) order(v, decreasing = FALSE)[2:(k + 1)]))
}

#' @keywords internal
.scprokar_build_knn_graph <- function(x, k = 15) {
  idx <- .scprokar_neighbor_index(x, k = k)
  if (ncol(idx) == 0) {
    return(Matrix::Matrix(0, nrow = nrow(x), ncol = nrow(x), sparse = TRUE))
  }
  i <- rep(seq_len(nrow(x)), each = ncol(idx))
  j <- as.vector(t(idx))
  adj <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(x), nrow(x)))
  adj <- (adj + Matrix::t(adj)) > 0
  methods::as(adj, "dgCMatrix")
}

#' @keywords internal
.scprokar_reduction_lookup <- function(sce, methods = NULL) {
  meta <- .scprokar_get_metadata(sce)
  integrations <- meta$integration$results
  mapping <- list()

  if (!is.null(integrations)) {
    for (nm in names(integrations)) {
      mapping[[nm]] <- integrations[[nm]]$reduction
    }
  }

  if (is.null(methods)) {
    if (length(mapping) > 0) {
      return(mapping)
    }
    names_out <- SingleCellExperiment::reducedDimNames(sce)
    out <- as.list(names_out)
    names(out) <- names_out
    return(out)
  }

  out <- list()
  for (item in methods) {
    if (!is.null(mapping[[item]])) {
      out[[item]] <- mapping[[item]]
    } else if (item %in% SingleCellExperiment::reducedDimNames(sce)) {
      out[[item]] <- item
    }
  }
  out
}

#' @keywords internal
.scprokar_adjusted_rand_index <- function(x, y) {
  x <- as.factor(x)
  y <- as.factor(y)
  tab <- table(x, y)
  n <- sum(tab)
  if (n <= 1) {
    return(NA_real_)
  }
  choose2 <- function(v) v * (v - 1) / 2
  sum_ij <- sum(choose2(tab))
  sum_i <- sum(choose2(rowSums(tab)))
  sum_j <- sum(choose2(colSums(tab)))
  expected <- sum_i * sum_j / choose2(n)
  max_index <- (sum_i + sum_j) / 2
  denom <- max_index - expected
  if (isTRUE(all.equal(denom, 0))) {
    return(NA_real_)
  }
  (sum_ij - expected) / denom
}

#' @keywords internal
.scprokar_entropy <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) {
    return(0)
  }
  -sum(p * log(p))
}

#' @keywords internal
.scprokar_nmi <- function(x, y) {
  x <- as.factor(x)
  y <- as.factor(y)
  joint <- prop.table(table(x, y))
  px <- rowSums(joint)
  py <- colSums(joint)
  mutual_info <- 0
  for (i in seq_len(nrow(joint))) {
    for (j in seq_len(ncol(joint))) {
      if (joint[i, j] > 0) {
        mutual_info <- mutual_info + joint[i, j] * log(joint[i, j] / (px[i] * py[j]))
      }
    }
  }
  hx <- .scprokar_entropy(px)
  hy <- .scprokar_entropy(py)
  if (isTRUE(all.equal(hx + hy, 0))) {
    return(NA_real_)
  }
  (2 * mutual_info) / (hx + hy)
}
