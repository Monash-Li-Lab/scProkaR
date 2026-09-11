#' Create a standardized bacterial single-cell object
#'
#' `CreateBacObject()` standardizes raw microbial single-cell counts into a
#' `SingleCellExperiment`, aligns cell and feature metadata, and records package
#' provenance used by downstream scProkaR workflows.
#'
#' @param x A gene-by-cell count matrix, a `SingleCellExperiment`, or a loaded
#'   Seurat object from an `.rds` file, or a Cell Ranger output directory.
#' @param counts_assay Assay name to use when `x` is a `SingleCellExperiment`.
#' @param cell_metadata Optional cell-level metadata with one row per cell.
#' @param feature_metadata Optional feature-level metadata with one row per
#'   gene.
#' @param feature_name_col Optional `rowData` column to use as the primary
#'   feature name when `x` is a `SingleCellExperiment` or 10x-derived object.
#'   If `NULL`, scProkaR will prefer common columns such as `Symbol`,
#'   `gene_name`, or `feature_name` when available. For 10x-derived input,
#'   `Symbol` is preferred by default when present.
#' @param sample_col Optional sample identifier column in `cell_metadata`.
#' @param batch_col Optional batch identifier column in `cell_metadata`.
#' @param condition_col Optional condition column in `cell_metadata`.
#' @param time_col Optional time or ordering column in `cell_metadata`.
#' @param sample_id_value Optional fixed sample identifier to assign to all
#'   cells during object creation. If `sample_col` is `NULL`, a `sample_id`
#'   column is created automatically.
#' @param batch_value Optional fixed batch identifier to assign to all cells
#'   during object creation. If `batch_col` is `NULL`, a `batch` column is
#'   created automatically.
#' @param condition_value Optional fixed condition value to assign to all cells
#'   during object creation. If `condition_col` is `NULL`, a `condition` column
#'   is created automatically.
#' @param time_value Optional fixed time value to assign to all cells during
#'   object creation. If `time_col` is `NULL`, a `time` column is created
#'   automatically.
#' @param organism Organism label stored in package metadata.
#' @param seurat_assay Assay name to extract when `x` is a Seurat object. If
#'   `NULL`, the active/default Seurat assay is used.
#' @param seurat_layer Layer to use as raw counts when `x` is a Seurat object.
#'   Defaults to `"counts"`.
#' @param transfer_reductions Logical indicating whether existing Seurat
#'   dimensional reductions (for example `pca` and `umap`) should be copied into
#'   `reducedDims(sce)`.
#' @param run_unintegrated Logical indicating whether to run a standard
#'   unintegrated preprocessing workflow after object creation. This computes
#'   normalized counts, identifies highly variable features, runs PCA, clusters
#'   cells, and computes an unintegrated UMAP without requiring Seurat.
#' @param unintegrated_feature_set Feature selection strategy passed to the
#'   unintegrated PCA step. Defaults to `"hvg"`.
#' @param unintegrated_dims Integer vector of dimensions to retain for PCA,
#'   neighbor graph construction, clustering, and UMAP.
#' @param unintegrated_cluster_col Column name used to store unintegrated
#'   cluster labels in `colData(sce)`.
#' @param unintegrated_umap_name Reduced-dimension name used to store the
#'   unintegrated UMAP embedding.
#' @param unintegrated_k Number of nearest neighbours used for unintegrated
#'   graph construction and clustering.
#' @param unintegrated_resolution Cluster granularity used for unintegrated
#'   clustering. Higher values typically produce more clusters.
#' @param unintegrated_algorithm Graph clustering backend for the unintegrated
#'   workflow.
#'
#'
#' For Cell Ranger / 10x input, users can pass the directory path directly.
#' scProkaR will import the matrix, standardize cell metadata, and by default
#' use `rowData(sce)$Symbol` as the feature name column when that column is
#' available. Set `feature_name_col` explicitly to override this behavior.
#'
#' @return A standardized `SingleCellExperiment`.
#' @export
#' @examples
#' set.seed(1)
#' genes <- paste0("gene", seq_len(40))
#' cells <- paste0("cell", seq_len(60))
#' counts <- matrix(
#'     rpois(length(genes) * length(cells), lambda = 3),
#'     nrow = length(genes),
#'     dimnames = list(genes, cells)
#' )
#' cell_meta <- data.frame(
#'     sample_id = rep(c("s1", "s2"), each = 30),
#'     time = rep(c(0, 30, 60), times = 20),
#'     row.names = cells
#' )
#' feature_meta <- data.frame(
#'     gene_class = rep(c("core", "accessory"), length.out = 40),
#'     row.names = genes
#' )
#'
#' sce <- CreateBacObject(
#'     counts,
#'     cell_metadata = cell_meta,
#'     feature_metadata = feature_meta,
#'     sample_col = "sample_id",
#'     time_col = "time",
#'     organism = "Escherichia coli"
#' )
#' sce
#' table(sce$sample_id, sce$time)
#' S4Vectors::metadata(sce)$scProkaR$organism
CreateBacObject <- function(
    x,
    counts_assay = "counts",
    cell_metadata = NULL,
    feature_metadata = NULL,
    feature_name_col = NULL,
    sample_col = NULL,
    batch_col = NULL,
    condition_col = NULL,
    time_col = NULL,
    sample_id_value = NULL,
    batch_value = NULL,
    condition_value = NULL,
    time_value = NULL,
    organism = "bacteria",
    seurat_assay = NULL,
    seurat_layer = "counts",
    transfer_reductions = TRUE,
    run_unintegrated = FALSE,
    unintegrated_feature_set = c("hvg", "all"),
    unintegrated_dims = seq_len(30),
    unintegrated_cluster_col = "unintegrated_clusters",
    unintegrated_umap_name = "umap.unintegrated",
    unintegrated_k = 20,
    unintegrated_resolution = 2,
    unintegrated_algorithm = c("louvain", "walktrap", "leiden")
) {
    unintegrated_feature_set <- match.arg(unintegrated_feature_set)
    unintegrated_algorithm <- match.arg(unintegrated_algorithm)

    ingested <- .scprokar_ingest_object_input(
        x,
        counts_assay = counts_assay,
        feature_name_col = feature_name_col,
        seurat_assay = seurat_assay,
        seurat_layer = seurat_layer,
        transfer_reductions = transfer_reductions
    )
    sce <- ingested$sce
    seurat_info <- ingested$seurat
    tenx_info <- ingested$tenx

    if (!is.null(sample_id_value) && is.null(sample_col))
        sample_col <- "sample_id"
    if (!is.null(batch_value) && is.null(batch_col)) batch_col <- "batch"
    if (!is.null(condition_value) && is.null(condition_col))
        condition_col <- "condition"
    if (!is.null(time_value) && is.null(time_col)) time_col <- "time"

    sce <- .scprokar_attach_standard_annotations(
        sce,
        cell_metadata = cell_metadata,
        feature_metadata = feature_metadata,
        sample_col = sample_col,
        batch_col = batch_col,
        condition_col = condition_col,
        time_col = time_col,
        sample_id_value = sample_id_value,
        batch_value = batch_value,
        condition_value = condition_value,
        time_value = time_value
    )

    sce <- .scprokar_record_creation_metadata(
        sce,
        organism = organism,
        sample_col = sample_col,
        batch_col = batch_col,
        condition_col = condition_col,
        time_col = time_col,
        seurat_info = seurat_info,
        tenx_info = tenx_info
    )

    if (isTRUE(run_unintegrated)) {
        sce <- .scprokar_apply_unintegrated_workflow(
            sce,
            feature_set = unintegrated_feature_set,
            dims = unintegrated_dims,
            cluster_col = unintegrated_cluster_col,
            umap_name = unintegrated_umap_name,
            k = unintegrated_k,
            resolution = unintegrated_resolution,
            algorithm = unintegrated_algorithm
        )
    }

    sce
}

#' @keywords internal
.scprokar_apply_fixed_metadata <- function(
    coldata,
    ids,
    sample_col = NULL,
    batch_col = NULL,
    condition_col = NULL,
    time_col = NULL,
    sample_id_value = NULL,
    batch_value = NULL,
    condition_value = NULL,
    time_value = NULL
) {
    out <- as.data.frame(coldata)
    if (is.null(rownames(out))) {
        rownames(out) <- ids
    }

    assignments <- list(
        list(column = sample_col, value = sample_id_value),
        list(column = batch_col, value = batch_value),
        list(column = condition_col, value = condition_value),
        list(column = time_col, value = time_value)
    )

    for (item in assignments) {
        if (is.null(item$column) || is.null(item$value)) {
            next
        }
        out[[item$column]] <- rep(item$value, length(ids))
    }

    out
}

#' Merge multiple scProkaR or SingleCellExperiment objects
#'
#' `MergeBacObjects()` combines multiple per-sample `SingleCellExperiment`
#' objects without relying on external `cbind()` method dispatch. This is useful
#' when merging sample-level objects created by [CreateBacObject()] before batch
#' integration.
#'
#' By default the function keeps the intersection of genes across inputs, which
#' is the safest choice before integration.
#'
#' @param ... Individual `SingleCellExperiment` objects to merge.
#' @param objects Optional list of `SingleCellExperiment` objects.
#' @param gene_mode Either `"intersect"` or `"union"`.
#'
#' @return A merged `SingleCellExperiment`.
#' @export
#' @examples
#' set.seed(1)
#' genes_a <- paste0("gene", 1:30)
#' genes_b <- paste0("gene", 11:40)
#' cells <- paste0("cell", seq_len(20))
#' counts_a <- matrix(
#'     rpois(length(genes_a) * 20, lambda = 4),
#'     nrow = length(genes_a),
#'     dimnames = list(genes_a, cells)
#' )
#' counts_b <- matrix(
#'     rpois(length(genes_b) * 20, lambda = 4),
#'     nrow = length(genes_b),
#'     dimnames = list(genes_b, cells)
#' )
#'
#' sce_a <- CreateBacObject(counts_a, sample_id_value = "sampleA")
#' sce_b <- CreateBacObject(counts_b, sample_id_value = "sampleB")
#'
#' ## Keep only the genes shared by both samples.
#' merged <- MergeBacObjects(sce_a, sce_b, gene_mode = "intersect")
#' dim(merged)
#' table(merged$sample_id)
#' head(colnames(merged), 3)
MergeBacObjects <- function(
    ...,
    objects = NULL,
    gene_mode = c("intersect", "union")
) {
    gene_mode <- match.arg(gene_mode)
    dots <- .scprokar_collect_merge_inputs(list(...), objects)
    dots <- .scprokar_prepare_objects_for_merge(dots)
    genes <- .scprokar_resolve_merge_genes(dots, gene_mode = gene_mode)

    merged <- .scprokar_assemble_merged_sce(dots, genes)
    merged <- .scprokar_add_merged_logcounts(merged, dots, genes)

    merged_reduction <- .scprokar_merge_reduced_dims2(dots)
    for (reduction_name in names(merged_reduction$kept)) {
        reduction_matrix <- merged_reduction$kept[[reduction_name]]
        SingleCellExperiment::reducedDim(merged, reduction_name) <-
            reduction_matrix
    }

    meta_list <- lapply(dots, .scprokar_get_metadata)
    merged_meta <- meta_list[[1]]
    merged_meta$merge <- list(
        n_objects = length(dots),
        gene_mode = gene_mode,
        n_genes = length(genes),
        n_cells = ncol(merged),
        cell_ids_unique = TRUE,
        reduced_dims = merged_reduction
    )
    .scprokar_set_metadata(merged, merged_meta)
}

#' @keywords internal
.scprokar_from_seurat <- function(
    x,
    seurat_assay = NULL,
    seurat_layer = "counts",
    transfer_reductions = TRUE
) {
    .scprokar_require("SeuratObject", "CreateBacObject() on Seurat input")

    if (is.null(seurat_assay)) {
        seurat_assay <- SeuratObject::DefaultAssay(x)
    }

    counts <- .scprokar_extract_seurat_layer(
        x,
        assay = seurat_assay,
        layer = seurat_layer
    )
    .scprokar_stopifnot_counts(counts)

    cell_metadata <- as.data.frame(x[[]])
    if (is.null(rownames(cell_metadata))) {
        rownames(cell_metadata) <- colnames(counts)
    }

    feature_metadata <- tryCatch(
        as.data.frame(x[[seurat_assay]][[]]),
        error = function(e) data.frame(row.names = rownames(counts))
    )
    feature_metadata <- .scprokar_align_data_frame(
        feature_metadata,
        ids = rownames(counts),
        what = "feature_metadata"
    )

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = .scprokar_as_dgC(counts)),
        colData = S4Vectors::DataFrame(.scprokar_align_data_frame(
            cell_metadata,
            colnames(counts),
            "cell_metadata"
        )),
        rowData = S4Vectors::DataFrame(feature_metadata)
    )

    normalized <- tryCatch(
        .scprokar_extract_seurat_layer(x, assay = seurat_assay, layer = "data"),
        error = function(e) NULL
    )
    if (!is.null(normalized) &&
        nrow(normalized) == nrow(sce) &&
        ncol(normalized) == ncol(sce)) {
        SummarizedExperiment::assay(sce, "logcounts") <- .scprokar_as_dgC(
            normalized[rownames(sce), colnames(sce), drop = FALSE]
        )
    }

    reduction_names <- character(0)
    if (isTRUE(transfer_reductions)) {
        reduction_names <- tryCatch(
            SeuratObject::Reductions(object = x),
            error = function(e) character(0)
        )
        for (reduction_name in reduction_names) {
            emb <- tryCatch(
                SeuratObject::Embeddings(
                    object = x,
                    reduction = reduction_name
                ),
                error = function(e) NULL
            )
            if (!is.null(emb)) {
                emb <- as.matrix(emb)
                if (all(colnames(sce) %in% rownames(emb))) {
                    emb <- emb[colnames(sce), , drop = FALSE]
                    SingleCellExperiment::reducedDim(sce, reduction_name) <- emb
                }
            }
        }
    }

    list(
        sce = sce,
        seurat = list(
            assay = seurat_assay,
            layer = seurat_layer,
            transferred_reductions = reduction_names
        )
    )
}

#' @keywords internal
.scprokar_extract_seurat_layer <- function(x, assay, layer) {
    data <- tryCatch(
        SeuratObject::LayerData(object = x, assay = assay, layer = layer),
        error = function(e) NULL
    )
    if (is.null(data)) {
        data <- tryCatch(
            SeuratObject::GetAssayData(
                object = x,
                assay = assay,
                layer = layer
            ),
            error = function(e) NULL
        )
    }
    if (is.null(data)) {
        data <- tryCatch(
            SeuratObject::GetAssayData(object = x, assay = assay, slot = layer),
            error = function(e) NULL
        )
    }
    if (is.null(data)) {
        stop(
            "Could not extract layer '", layer, "' from Seurat assay '", assay,
            "'. Check the assay/layer names in the loaded object.",
            call. = FALSE
        )
    }
    data
}

#' @keywords internal
.scprokar_from_10x_dir <- function(x) {
    source_dir <- normalizePath(x, winslash = "/", mustWork = TRUE)

    if (requireNamespace("DropletUtils", quietly = TRUE)) {
        sce <- tryCatch(
            DropletUtils::read10xCounts(source_dir),
            error = function(e) NULL
        )
        if (!is.null(sce)) {
            sce <- .scprokar_ensure_cell_names(sce)
            counts <- SummarizedExperiment::assay(sce, "counts")
            SummarizedExperiment::assay(sce, "counts") <-
                .scprokar_as_dgC(counts)
            return(list(
                sce = sce,
                tenx = list(
                    source_dir = source_dir,
                    reader = "DropletUtils::read10xCounts"
                )
            ))
        }
    }

    matrix_dir <- .scprokar_locate_10x_matrix_dir(source_dir)
    counts <- .scprokar_read_10x_matrix(matrix_dir)
    barcodes <- .scprokar_read_10x_table(matrix_dir, "barcodes")
    features <- .scprokar_read_10x_table(matrix_dir, "features")

    barcodes <- as.character(barcodes[[1]])
    feature_ids <- as.character(features[[1]])
    feature_names <- if (ncol(features) >= 2) {
        as.character(features[[2]])
    } else {
        feature_ids
    }
    rownames(counts) <- make.unique(feature_names)
    colnames(counts) <- barcodes

    feature_df <- data.frame(
        feature_id = feature_ids,
        feature_name = feature_names,
        stringsAsFactors = FALSE,
        row.names = rownames(counts)
    )
    if (ncol(features) >= 3) {
        feature_df$feature_type <- as.character(features[[3]])
    }

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = .scprokar_as_dgC(counts)),
        rowData = S4Vectors::DataFrame(feature_df),
        colData = S4Vectors::DataFrame(row.names = barcodes)
    )

    list(
        sce = sce,
        tenx = list(
            source_dir = source_dir,
            matrix_dir = matrix_dir,
            reader = "manual_matrix_market"
        )
    )
}

#' @keywords internal
.scprokar_ensure_cell_names <- function(sce) {
    if (!is.null(colnames(sce)) && all(nzchar(colnames(sce)))) {
        return(sce)
    }

    cd <- as.data.frame(SummarizedExperiment::colData(sce))
    barcode_col <- intersect(
        c("Barcode", "barcode", "cell", "cell_id"),
        colnames(cd)
    )

    if (length(barcode_col) >= 1L) {
        candidate <- as.character(cd[[barcode_col[[1]]]])
        if (length(candidate) == ncol(sce) && all(nzchar(candidate))) {
            colnames(sce) <- make.unique(candidate)
            return(sce)
        }
    }

    if (!is.null(rownames(cd)) &&
        length(rownames(cd)) == ncol(sce) &&
        all(nzchar(rownames(cd)))) {
        colnames(sce) <- make.unique(rownames(cd))
        return(sce)
    }

    stop(
        "Could not determine cell barcodes from the imported 10x object. ",
        "Expected barcodes in `colnames(sce)` or in a colData column such as ",
        "`Barcode`.",
        call. = FALSE
    )
}

#' @keywords internal
.scprokar_ensure_feature_names <- function(sce, counts_assay = "counts") {
    current_names <- rownames(sce)
    if (!is.null(current_names) && all(nzchar(current_names))) {
        return(sce)
    }

    rd <- as.data.frame(SummarizedExperiment::rowData(sce))
    candidate_cols <- intersect(
        c(
            "Symbol", "symbol", "gene_name", "feature_name",
            "ID", "id", "gene_id"
        ),
        colnames(rd)
    )

    for (candidate_col in candidate_cols) {
        candidate <- as.character(rd[[candidate_col]])
        if (length(candidate) == nrow(sce) && any(nzchar(candidate))) {
            rownames(sce) <- make.unique(ifelse(
                is.na(candidate) | !nzchar(candidate),
                paste0("feature_", seq_along(candidate)),
                candidate
            ))
            assay_names <- SummarizedExperiment::assayNames(sce)
            for (assay_name in assay_names) {
                assay_mat <- SummarizedExperiment::assay(sce, assay_name)
                rownames(assay_mat) <- rownames(sce)
                SummarizedExperiment::assay(
                    sce, assay_name,
                    withDimnames = FALSE
                ) <- assay_mat
            }
            return(sce)
        }
    }

    counts <- SummarizedExperiment::assay(sce, counts_assay)
    if (!is.null(rownames(counts)) && all(nzchar(rownames(counts)))) {
        rownames(sce) <- make.unique(rownames(counts))
        return(sce)
    }

    stop(
        "Could not determine feature names from the imported object. ",
        "Expected names in `rownames(sce)` or in rowData columns such as ",
        "`Symbol` or `ID`.",
        call. = FALSE
    )
}

#' @keywords internal
.scprokar_promote_feature_names <- function(sce, feature_name_col = NULL) {
    rd <- as.data.frame(SummarizedExperiment::rowData(sce))
    candidate_cols <- if (!is.null(feature_name_col)) {
        feature_name_col
    } else {
        intersect(
            c(
                "Symbol", "symbol", "gene_name", "feature_name",
                "ID", "id", "gene_id"
            ),
            colnames(rd)
        )
    }

    for (candidate_col in candidate_cols) {
        if (!candidate_col %in% colnames(rd)) {
            next
        }
        candidate <- as.character(rd[[candidate_col]])
        if (length(candidate) != nrow(sce) || !any(nzchar(candidate))) {
            next
        }
        replacement <- ifelse(
            is.na(candidate) | !nzchar(candidate),
            if (!is.null(rownames(sce)) &&
                length(rownames(sce)) == nrow(sce)) rownames(sce) else
                paste0("feature_", seq_len(nrow(sce))),
            candidate
        )
        replacement <- make.unique(replacement)

        current_names <- rownames(sce)
        if (!is.null(current_names) && identical(current_names, replacement)) {
            return(sce)
        }

        if (!"feature_id" %in% colnames(rd) && !is.null(current_names)) {
            rd$feature_id <- current_names
        }
        rownames(rd) <- replacement
        SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(rd)
        rownames(sce) <- replacement

        assay_names <- SummarizedExperiment::assayNames(sce)
        for (assay_name in assay_names) {
            assay_mat <- SummarizedExperiment::assay(sce, assay_name)
            rownames(assay_mat) <- replacement
            SummarizedExperiment::assay(
                sce, assay_name,
                withDimnames = FALSE
            ) <- assay_mat
        }
        return(sce)
    }

    sce
}

#' @keywords internal
.scprokar_prepare_sce_input <- function(
    sce,
    counts_assay = "counts",
    feature_name_col = NULL
) {
    if (!counts_assay %in% SummarizedExperiment::assayNames(sce)) {
        stop("Assay '", counts_assay, "' was not found in `x`.", call. = FALSE)
    }

    sce <- .scprokar_ensure_cell_names(sce)
    sce <- .scprokar_promote_feature_names(
        sce,
        feature_name_col = feature_name_col
    )
    sce <- .scprokar_ensure_feature_names(sce, counts_assay = counts_assay)

    assay_names <- SummarizedExperiment::assayNames(sce)
    for (assay_name in assay_names) {
        assay_mat <- SummarizedExperiment::assay(sce, assay_name)
        if (is.null(colnames(assay_mat)) ||
            !identical(colnames(assay_mat), colnames(sce))) {
            colnames(assay_mat) <- colnames(sce)
        }
        if (is.null(rownames(assay_mat)) ||
            !identical(rownames(assay_mat), rownames(sce))) {
            rownames(assay_mat) <- rownames(sce)
        }
        SummarizedExperiment::assay(
            sce, assay_name,
            withDimnames = FALSE
        ) <- assay_mat
    }

    sce
}

#' @keywords internal
.scprokar_locate_10x_matrix_dir <- function(source_dir) {
    candidates <- c(
        source_dir,
        file.path(source_dir, "filtered_feature_bc_matrix"),
        file.path(source_dir, "raw_feature_bc_matrix"),
        file.path(source_dir, "outs", "filtered_feature_bc_matrix"),
        file.path(source_dir, "outs", "raw_feature_bc_matrix")
    )
    candidates <- unique(candidates[dir.exists(candidates)])

    for (candidate in candidates) {
        has_matrix <- any(file.exists(file.path(
            candidate,
            c("matrix.mtx", "matrix.mtx.gz")
        )))
        has_barcodes <- any(file.exists(file.path(
            candidate,
            c("barcodes.tsv", "barcodes.tsv.gz")
        )))
        has_features <- any(file.exists(file.path(
            candidate,
            c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")
        )))
        if (has_matrix && has_barcodes && has_features) {
            return(candidate)
        }
    }

    stop(
        "Could not locate a 10x matrix directory under '", source_dir,
        "'. Expected files like matrix.mtx(.gz), barcodes.tsv(.gz), and ",
        "features.tsv(.gz).",
        call. = FALSE
    )
}

#' @keywords internal
.scprokar_read_10x_matrix <- function(matrix_dir) {
    matrix_file <- .scprokar_first_existing_file(
        matrix_dir,
        c("matrix.mtx", "matrix.mtx.gz")
    )
    con <- if (grepl("\\.gz$", matrix_file))
        gzfile(matrix_file, open = "rt") else file(matrix_file, open = "rt")
    on.exit(close(con), add = TRUE)
    Matrix::readMM(con)
}

#' @keywords internal
.scprokar_read_10x_table <- function(
    matrix_dir,
    kind = c("barcodes", "features")
) {
    kind <- match.arg(kind)
    choices <- switch(kind,
        barcodes = c("barcodes.tsv", "barcodes.tsv.gz"),
        features = c(
            "features.tsv", "features.tsv.gz",
            "genes.tsv", "genes.tsv.gz"
        )
    )
    table_file <- .scprokar_first_existing_file(matrix_dir, choices)
    con <- if (grepl("\\.gz$", table_file))
        gzfile(table_file, open = "rt") else file(table_file, open = "rt")
    on.exit(close(con), add = TRUE)
    utils::read.delim(con, header = FALSE, stringsAsFactors = FALSE)
}

#' @keywords internal
.scprokar_first_existing_file <- function(path, candidates) {
    full_paths <- file.path(path, candidates)
    hits <- full_paths[file.exists(full_paths)]
    if (!length(hits)) {
        stop(
            "None of the expected files were found in '", path, "': ",
            paste(candidates, collapse = ", "),
            call. = FALSE
        )
    }
    hits[[1]]
}

#' @keywords internal
.scprokar_expand_counts <- function(mat, genes) {
    out <- Matrix::Matrix(
        0,
        nrow = length(genes),
        ncol = ncol(mat),
        sparse = TRUE
    )
    rownames(out) <- genes
    colnames(out) <- colnames(mat)
    idx <- match(rownames(mat), genes)
    keep <- !is.na(idx)
    if (any(keep)) {
        out[idx[keep], ] <- mat[keep, , drop = FALSE]
    }
    .scprokar_as_dgC(out)
}

#' @keywords internal
.scprokar_merge_rowdata <- function(rowdata_list, genes) {
    all_columns <- unique(unlist(
        lapply(rowdata_list, colnames),
        use.names = FALSE
    ))
    out <- data.frame(row.names = genes)
    for (col in all_columns) {
        values <- rep(NA, length(genes))
        for (df in rowdata_list) {
            if (!col %in% colnames(df)) {
                next
            }
            idx <- match(rownames(df), genes)
            fill <- !is.na(idx) & is.na(values[idx])
            if (any(fill)) {
                values[idx[fill]] <- df[[col]][fill]
            }
        }
        out[[col]] <- values
    }
    out
}

#' @keywords internal
.scprokar_merge_reduced_dims2 <- function(objects) {
    reduction_names <- lapply(objects, SingleCellExperiment::reducedDimNames)
    common_names <- Reduce(intersect, reduction_names)
    if (length(common_names) == 0L) {
        return(list(
            kept = list(),
            skipped = list(),
            n_objects = length(objects)
        ))
    }

    kept <- list()
    skipped <- list()
    cell_ids <- unlist(lapply(objects, colnames), use.names = FALSE)

    for (reduction_name in common_names) {
        reductions <- lapply(objects, function(sce) {
            SingleCellExperiment::reducedDim(sce, reduction_name)
        })
        ndim <- vapply(reductions, ncol, integer(1))
        if (length(unique(ndim)) != 1L) {
            skipped[[reduction_name]] <- sprintf(
                "dimension mismatch (ncol: %s)",
                paste(ndim, collapse = ",")
            )
            next
        }

        aligned <- lapply(seq_along(objects), function(i) {
            reduction <- as.matrix(reductions[[i]])
            reduction_cells <- rownames(reduction)
            target_cells <- colnames(objects[[i]])

            if (is.null(reduction_cells) || anyNA(reduction_cells)) {
                if (nrow(reduction) == length(target_cells)) {
                    rownames(reduction) <- target_cells
                    return(reduction)
                }
                return(NULL)
            }

            if (all(target_cells %in% reduction_cells)) {
                return(reduction[target_cells, , drop = FALSE])
            }

            idx <- match(target_cells, reduction_cells)
            aligned_reduction <- matrix(
                NA_real_,
                nrow = length(target_cells),
                ncol = ncol(reduction),
                dimnames = list(target_cells, colnames(reduction))
            )
            hit <- !is.na(idx)
            if (!any(hit)) {
                return(NULL)
            }
            aligned_reduction[hit, ] <- reduction[idx[hit], , drop = FALSE]
            aligned_reduction
        })
        if (any(vapply(aligned, is.null, logical(1)))) {
            skipped[[reduction_name]] <-
                "unable to align reduced dimensions by cell IDs"
            next
        }

        merged_matrix <- do.call(rbind, aligned)
        rownames(merged_matrix) <- cell_ids
        kept[[reduction_name]] <- merged_matrix
    }

    list(kept = kept, skipped = skipped, n_objects = length(objects))
}

#' @keywords internal
.scprokar_prepare_objects_for_merge <- function(objects) {
    all_ids <- unlist(lapply(objects, colnames), use.names = FALSE)
    if (!anyDuplicated(all_ids)) {
        return(objects)
    }

    lapply(seq_along(objects), function(i) {
        sce <- objects[[i]]
        new_ids <- .scprokar_merge_cell_ids(sce, object_index = i)

        colnames(sce) <- new_ids
        cd <- as.data.frame(SummarizedExperiment::colData(sce))
        cd$original_cell_id <- rownames(cd)
        rownames(cd) <- new_ids
        SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cd)

        for (reduction_name in SingleCellExperiment::reducedDimNames(sce)) {
            emb <- SingleCellExperiment::reducedDim(sce, reduction_name)
            rownames(emb) <- new_ids
            SingleCellExperiment::reducedDim(sce, reduction_name) <- emb
        }

        meta <- .scprokar_get_metadata(sce)
        meta$merge <- c(meta$merge, list(renamed_for_merge = TRUE))
        .scprokar_set_metadata(sce, meta)
    })
}

#' @keywords internal
.scprokar_merge_cell_ids <- function(sce, object_index) {
    meta <- .scprokar_get_metadata(sce)
    cols <- meta$columns
    cd <- as.data.frame(SummarizedExperiment::colData(sce))

    prefix <- NULL
    for (candidate in c(
        cols$sample_col, cols$batch_col, "sample_id", "batch"
    )) {
        if (is.null(candidate) || !candidate %in% colnames(cd)) {
            next
        }
        values <- unique(as.character(cd[[candidate]]))
        values <- values[!is.na(values)]
        if (length(values) == 1) {
            prefix <- values[[1]]
            break
        }
    }

    if (is.null(prefix) || !nzchar(prefix)) {
        prefix <- paste0("sample", object_index)
    }

    proposed <- paste(prefix, colnames(sce), sep = "_")
    make.unique(proposed, sep = "_")
}

#' Import raw input into a standardized SingleCellExperiment
#'
#' Dispatches on the class of `x` and returns the imported object together
#' with any Seurat or 10x provenance recorded during import.
#'
#' @keywords internal
#' @noRd
.scprokar_ingest_object_input <- function(
    x,
    counts_assay = "counts",
    feature_name_col = NULL,
    seurat_assay = NULL,
    seurat_layer = "counts",
    transfer_reductions = TRUE
) {
    seurat_info <- NULL
    tenx_info <- NULL

    if (methods::is(x, "SingleCellExperiment")) {
        sce <- .scprokar_prepare_sce_input(
            x,
            counts_assay = counts_assay,
            feature_name_col = feature_name_col
        )
        if (!counts_assay %in% SummarizedExperiment::assayNames(sce)) {
            stop(
                "Assay '", counts_assay, "' was not found in `x`.",
                call. = FALSE
            )
        }
        counts <- SummarizedExperiment::assay(sce, counts_assay)
        .scprokar_stopifnot_counts(counts)
        SummarizedExperiment::assay(sce, "counts") <- .scprokar_as_dgC(counts)
    } else if (methods::is(x, "Seurat")) {
        seurat_payload <- .scprokar_from_seurat(
            x,
            seurat_assay = seurat_assay,
            seurat_layer = seurat_layer,
            transfer_reductions = transfer_reductions
        )
        sce <- seurat_payload$sce
        seurat_info <- seurat_payload$seurat
    } else if (is.character(x) && length(x) == 1L && dir.exists(x)) {
        tenx_payload <- .scprokar_from_10x_dir(x)
        sce <- .scprokar_prepare_sce_input(
            tenx_payload$sce,
            counts_assay = "counts",
            feature_name_col = feature_name_col
        )
        tenx_info <- tenx_payload$tenx
    } else {
        .scprokar_stopifnot_counts(x)
        counts <- .scprokar_as_dgC(x)
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = counts)
        )
    }

    list(sce = sce, seurat = seurat_info, tenx = tenx_info)
}

#' Align cell and feature metadata onto a SingleCellExperiment
#'
#' Aligns supplied (or existing) metadata to the object dimnames, applies any
#' fixed annotation values, and validates the requested annotation columns.
#'
#' @keywords internal
#' @noRd
.scprokar_attach_standard_annotations <- function(
    sce,
    cell_metadata = NULL,
    feature_metadata = NULL,
    sample_col = NULL,
    batch_col = NULL,
    condition_col = NULL,
    time_col = NULL,
    sample_id_value = NULL,
    batch_value = NULL,
    condition_value = NULL,
    time_value = NULL
) {
    cells <- colnames(sce)
    genes <- rownames(sce)

    merged_coldata <- .scprokar_align_data_frame(
        if (!is.null(cell_metadata)) cell_metadata else
            as.data.frame(SummarizedExperiment::colData(sce)),
        ids = cells,
        what = "cell_metadata"
    )
    merged_coldata <- .scprokar_apply_fixed_metadata(
        merged_coldata,
        ids = cells,
        sample_col = sample_col,
        batch_col = batch_col,
        condition_col = condition_col,
        time_col = time_col,
        sample_id_value = sample_id_value,
        batch_value = batch_value,
        condition_value = condition_value,
        time_value = time_value
    )
    merged_rowdata <- .scprokar_align_data_frame(
        if (!is.null(feature_metadata)) feature_metadata else
            as.data.frame(SummarizedExperiment::rowData(sce)),
        ids = genes,
        what = "feature_metadata"
    )

    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(merged_coldata)
    SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(merged_rowdata)

    .scprokar_match_columns(
        sce,
        columns = c(sample_col, batch_col, condition_col, time_col),
        label = "annotation"
    )

    sce
}

#' Record scProkaR provenance metadata on a newly created object
#'
#' @keywords internal
#' @noRd
.scprokar_record_creation_metadata <- function(
    sce,
    organism,
    sample_col = NULL,
    batch_col = NULL,
    condition_col = NULL,
    time_col = NULL,
    seurat_info = NULL,
    tenx_info = NULL
) {
    meta <- .scprokar_get_metadata(sce)
    meta$package <- "scProkaR"
    meta$version <- .scprokar_package_version()
    meta$created <- as.character(Sys.time())
    meta$organism <- organism
    meta$columns <- list(
        sample_col = sample_col,
        batch_col = batch_col,
        condition_col = condition_col,
        time_col = time_col
    )
    meta$assays <- list(counts = "counts")
    if (!is.null(seurat_info)) {
        meta$seurat <- seurat_info
    }
    if (!is.null(tenx_info)) {
        meta$tenx <- tenx_info
    }

    .scprokar_set_metadata(sce, meta)
}

#' Collect and validate the objects supplied to MergeBacObjects()
#'
#' @keywords internal
#' @noRd
.scprokar_collect_merge_inputs <- function(dots, objects = NULL) {
    if (!is.null(objects)) {
        dots <- c(dots, objects)
    }
    dots <- Filter(Negate(is.null), dots)

    if (length(dots) < 2) {
        stop(
            "Provide at least two SingleCellExperiment objects to merge.",
            call. = FALSE
        )
    }
    if (!all(vapply(dots, methods::is, logical(1), "SingleCellExperiment"))) {
        stop(
        "All inputs to MergeBacObjects() must be SingleCellExperiment objects.",
        call. = FALSE
        )
    }

    dots
}

#' Resolve the gene universe shared by the objects being merged
#'
#' @keywords internal
#' @noRd
.scprokar_resolve_merge_genes <- function(objects, gene_mode) {
    gene_sets <- lapply(objects, rownames)
    genes <- if (gene_mode == "intersect") {
        Reduce(intersect, gene_sets)
    } else {
        Reduce(union, gene_sets)
    }
    if (length(genes) == 0) {
        stop(
            "No genes remained after applying `gene_mode = \"", gene_mode,
            "\"`.",
            call. = FALSE
        )
    }
    genes
}

#' Build the merged SingleCellExperiment from counts and aligned metadata
#'
#' @keywords internal
#' @noRd
.scprokar_assemble_merged_sce <- function(objects, genes) {
    counts_list <- lapply(objects, function(sce) {
        counts <- SummarizedExperiment::assay(sce, "counts")
        .scprokar_expand_counts(counts, genes)
    })
    merged_counts <- Reduce(Matrix::cbind2, counts_list)
    merged_counts <- .scprokar_as_dgC(merged_counts)

    coldata_list <- lapply(objects, function(sce) {
        as.data.frame(SummarizedExperiment::colData(sce))
    })
    merged_coldata <- do.call(rbind, coldata_list)
    rownames(merged_coldata) <- unlist(
        lapply(objects, colnames),
        use.names = FALSE
    )

    rowdata_list <- lapply(objects, function(sce) {
        df <- as.data.frame(SummarizedExperiment::rowData(sce))
        df <- df[match(rownames(sce), rownames(df)), , drop = FALSE]
        rownames(df) <- rownames(sce)
        df
    })
    merged_rowdata <- .scprokar_merge_rowdata(rowdata_list, genes)

    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = merged_counts),
        colData = S4Vectors::DataFrame(merged_coldata),
        rowData = S4Vectors::DataFrame(merged_rowdata)
    )
}

#' Add a merged logcounts assay when every input object provides one
#'
#' @keywords internal
#' @noRd
.scprokar_add_merged_logcounts <- function(merged, objects, genes) {
    has_logcounts <- all(vapply(
        objects,
        function(sce) "logcounts" %in%
            SummarizedExperiment::assayNames(sce),
        logical(1)
    ))
    if (!has_logcounts) {
        return(merged)
    }

    logcounts_list <- lapply(objects, function(sce) {
        logcounts <- SummarizedExperiment::assay(sce, "logcounts")
        .scprokar_expand_counts(logcounts, genes)
    })
    SummarizedExperiment::assay(merged, "logcounts") <- .scprokar_as_dgC(
        Reduce(Matrix::cbind2, logcounts_list)
    )
    merged
}

#' Run the unintegrated workflow and record its parameters in metadata
#'
#' @keywords internal
#' @noRd
.scprokar_apply_unintegrated_workflow <- function(
    sce,
    feature_set,
    dims,
    cluster_col,
    umap_name,
    k,
    resolution,
    algorithm
) {
    sce <- .scprokar_run_unintegrated_workflow(
        sce,
        feature_set = feature_set,
        dims = dims,
        cluster_col = cluster_col,
        umap_name = umap_name,
        k = k,
        resolution = resolution,
        algorithm = algorithm
    )
    meta <- .scprokar_get_metadata(sce)
    meta$unintegrated <- list(
        feature_set = feature_set,
        dims = dims,
        cluster_col = cluster_col,
        umap_name = umap_name,
        k = k,
        resolution = resolution,
        algorithm = algorithm
    )
    .scprokar_set_metadata(sce, meta)
}
