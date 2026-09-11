#' Internal utilities for scProkaR
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
        ## stop() concatenates its arguments, so build the message from
        ## separate pieces rather than with paste(), which BiocCheck flags
        ## in condition signals.
        purpose <- if (is.null(reason)) "" else paste0(" for ", reason)
        stop(
            "Package '", pkg, "' is required", purpose,
            ". Install it first.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' @keywords internal
.scprokar_package_version <- function() {
    desc <- tryCatch(
        utils::packageDescription("scProkaR"),
        warning = function(w) NULL,
        error = function(e) NULL
    )
    if (!is.null(desc) && !is.null(desc$Version)) {
        return(as.character(desc$Version))
    }

    desc_file <- system.file("DESCRIPTION", package = "scProkaR")
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
            "`x` must be a matrix-like gene-by-cell count object or a ",
            "SingleCellExperiment.",
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
    if (is.matrix(x)) {
        x <- Matrix::Matrix(x, sparse = TRUE)
    }
    ## `is.logical()` is FALSE for every S4 Matrix, so test the class instead.
    ## Multiplying by 1 promotes a logical Matrix to a numeric one and avoids
    ## the deprecated direct coercion from lMatrix to dgCMatrix.
    if (inherits(x, "lMatrix")) {
        x <- x * 1
    }
    if (inherits(x, "Matrix")) {
        x <- methods::as(x, "CsparseMatrix")
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
        stop(
            "`", what, "` is missing identifiers: ",
            paste(missing_ids, collapse = ", "),
            call. = FALSE
        )
    }
    x[ids, , drop = FALSE]
}

#' @keywords internal
.scprokar_get_metadata <- function(sce) {
    meta <- S4Vectors::metadata(sce)$scProkaR
    if (is.null(meta)) {
        meta <- list()
    }
    meta
}

#' @keywords internal
.scprokar_set_metadata <- function(sce, meta) {
    all_meta <- S4Vectors::metadata(sce)
    all_meta$scProkaR <- meta
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
.scprokar_normalize_logcounts <- function(
    sce,
    assay_name = "counts",
    scale_factor = 1e4
) {
    counts <- SummarizedExperiment::assay(sce, assay_name)
    libsize <- Matrix::colSums(counts)
    scaling <- scale_factor / pmax(libsize, 1)
    norm <- counts %*% Matrix::Diagonal(x = scaling)
    lognorm <- .scprokar_as_dgC(log1p(norm))
    dimnames(lognorm) <- dimnames(counts)
    SummarizedExperiment::assay(sce, "logcounts", withDimnames = FALSE) <-
        lognorm
    sce
}

#' @keywords internal
.scprokar_select_features <- function(
    sce,
    feature_set = c("hvg", "all"),
    nfeatures = 2000
) {
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
        n_cells <- ncol(logx)
        center <- Matrix::rowMeans(logx)
        mean_sq <- Matrix::rowSums(logx^2) / max(n_cells, 1)
        variance <- pmax(as.numeric(mean_sq) - as.numeric(center)^2, 0)
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
    mat <- SummarizedExperiment::assay(
        sce, assay_name
    )[features, , drop = FALSE]
    x_sparse <- Matrix::t(.scprokar_as_dgC(mat))
    rank_k <- min(ncomponents, max(1, ncol(x_sparse) - 1))
    full_rank <- min(nrow(x_sparse), ncol(x_sparse))

    use_irlba <- requireNamespace("irlba", quietly = TRUE) &&
        ncol(x_sparse) > 2 &&
        rank_k < full_rank

    if (use_irlba) {
        center <- Matrix::colMeans(x_sparse)
        mean_sq <- Matrix::colMeans(x_sparse^2)
        scale_vec <- sqrt(pmax(mean_sq - center^2, 0))
        scale_vec[!is.finite(scale_vec) | scale_vec <= 0] <- 1

        pcs <- irlba::prcomp_irlba(
            x_sparse,
            n = rank_k,
            center = center,
            scale. = scale_vec
        )
    } else {
        if ((nrow(x_sparse) * ncol(x_sparse)) > 5e7) {
            stop(
                "Sparse-aware PCA requires the optional package 'irlba' for ",
                "larger datasets. ",
                "Install it or reduce the feature set before running ",
                "`run_unintegrated = TRUE`.",
                call. = FALSE
            )
        }
        x_dense <- as.matrix(x_sparse)
        x_dense <- scale(x_dense, center = TRUE, scale = TRUE)
        x_dense[is.na(x_dense)] <- 0
        pcs <- stats::prcomp(
            x_dense,
            rank. = rank_k,
            center = FALSE,
            scale. = FALSE
        )
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
    dims = seq_len(30),
    reduction_name = "PCA"
) {
    feature_set <- match.arg(feature_set)
    if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
        sce <- .scprokar_normalize_logcounts(sce, assay_name = assay_name)
    }
    if (!reduction_name %in% SingleCellExperiment::reducedDimNames(sce) ||
        max(dims) > ncol(
            SingleCellExperiment::reducedDim(sce, reduction_name)
        )) {
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
.scprokar_run_unintegrated_workflow <- function(
    sce,
    feature_set = "hvg",
    dims = seq_len(30),
    cluster_col = "unintegrated_clusters",
    umap_name = "umap.unintegrated",
    k = 20,
    resolution = 2,
    algorithm = "louvain"
) {
    sce <- .scprokar_prepare_reduction(
        sce,
        assay_name = "counts",
        feature_set = feature_set,
        dims = dims,
        reduction_name = "PCA"
    )
    sce <- RunIntegratedClustering(
        sce,
        reduction = "PCA",
        cluster_col = cluster_col,
        dims = dims,
        k = k,
        algorithm = algorithm,
        resolution = resolution
    )

    if (requireNamespace("uwot", quietly = TRUE)) {
        sce <- RunIntegratedUMAP(
            sce,
            reduction = "PCA",
            umap_name = umap_name,
            n_neighbors = k
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
    adj <- adj * 1
    methods::as(adj, "dgCMatrix")
}

#' @keywords internal
.scprokar_reduction_lookup <- function(sce, methods = NULL) {
    meta <- .scprokar_get_metadata(sce)
    integrations <- meta$integration$results
    mapping <- list()
    available <- SingleCellExperiment::reducedDimNames(sce)

    if (!is.null(integrations)) {
        for (nm in names(integrations)) {
            mapping[[nm]] <- integrations[[nm]]$reduction
        }
    }

    ## Aliases inferred from an `integrated_*` name that was never registered.
    ## They are a convenience, so they rank BELOW an exact reduced-dimension
    ## name: asking for a reduction by its own name must never silently score
    ## a different one.
    discovered <- list()
    registered <- mapping
    registered_reductions <- unname(unlist(registered, use.names = FALSE))
    for (reduction_name in grep("^integrated_", available, value = TRUE)) {
        if (reduction_name %in% registered_reductions) {
            next
        }
        method_name <- .scprokar_infer_integration_method_name(reduction_name)
        if (!nzchar(method_name) || !is.null(registered[[method_name]]) ||
                !is.null(discovered[[method_name]])) {
            method_name <- reduction_name
        }
        discovered[[method_name]] <- reduction_name
    }
    mapping <- c(
        registered,
        discovered[setdiff(names(discovered), names(registered))]
    )

    if (is.null(methods)) {
        if (length(mapping) > 0) {
            return(mapping)
        }
        names_out <- available
        out <- as.list(names_out)
        names(out) <- names_out
        return(out)
    }

    out <- list()
    for (item in methods) {
        if (!is.null(registered[[item]])) {
            out[[item]] <- registered[[item]]
        } else if (item %in% available) {
            out[[item]] <- item
        } else if (!is.null(discovered[[item]])) {
            out[[item]] <- discovered[[item]]
        } else if (tolower(item) %in% tolower(names(mapping))) {
            matched <- names(mapping)[match(
                tolower(item), tolower(names(mapping))
            )]
            out[[matched]] <- mapping[[matched]]
        } else if (tolower(item) %in% tolower(available)) {
            matched <- available[match(tolower(item), tolower(available))]
            out[[matched]] <- matched
        }
    }
    out
}

#' @keywords internal
.scprokar_adjusted_rand_index <- function(x, y) {
    x <- as.factor(x)
    y <- as.factor(y)
    if (length(x) <= 1L) {
        return(NA_real_)
    }
    ## Delegate to bluster rather than carrying our own implementation of a
    ## standard statistic. Verified to agree to machine precision with the
    ## previous hand-written version.
    ari <- bluster::pairwiseRand(x, y, mode = "index", adjusted = TRUE)
    if (!is.finite(ari)) {
        return(NA_real_)
    }
    as.numeric(ari)
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
                mutual_info <- mutual_info +
                    joint[i, j] * log(joint[i, j] / (px[i] * py[j]))
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
