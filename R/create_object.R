#' Create a standardized bacterial single-cell object
#'
#' `CreateBacObject()` standardizes raw microbial single-cell counts into a
#' `SingleCellExperiment`, aligns cell and feature metadata, and records package
#' provenance used by downstream SCProkaR workflows.
#'
#' @param x A gene-by-cell count matrix, a `SingleCellExperiment`, or a loaded
#'   Seurat object from an `.rds` file.
#' @param counts_assay Assay name to use when `x` is a `SingleCellExperiment`.
#' @param cell_metadata Optional cell-level metadata with one row per cell.
#' @param feature_metadata Optional feature-level metadata with one row per gene.
#' @param sample_col Optional sample identifier column in `cell_metadata`.
#' @param batch_col Optional batch identifier column in `cell_metadata`.
#' @param condition_col Optional condition column in `cell_metadata`.
#' @param time_col Optional time or ordering column in `cell_metadata`.
#' @param organism Organism label stored in package metadata.
#' @param seurat_assay Assay name to extract when `x` is a Seurat object. If
#'   `NULL`, the active/default Seurat assay is used.
#' @param seurat_layer Layer to use as raw counts when `x` is a Seurat object.
#'   Defaults to `"counts"`.
#' @param transfer_reductions Logical indicating whether existing Seurat
#'   dimensional reductions (for example `pca` and `umap`) should be copied into
#'   `reducedDims(sce)`.
#'
#' @return A standardized `SingleCellExperiment`.
#' @export
CreateBacObject <- function(
    x,
    counts_assay = "counts",
    cell_metadata = NULL,
    feature_metadata = NULL,
    sample_col = NULL,
    batch_col = NULL,
    condition_col = NULL,
    time_col = NULL,
    organism = "bacteria",
    seurat_assay = NULL,
    seurat_layer = "counts",
    transfer_reductions = TRUE
) {
  seurat_info <- NULL

  if (methods::is(x, "SingleCellExperiment")) {
    sce <- x
    if (!counts_assay %in% SummarizedExperiment::assayNames(sce)) {
      stop("Assay '", counts_assay, "' was not found in `x`.", call. = FALSE)
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
  } else {
    .scprokar_stopifnot_counts(x)
    counts <- .scprokar_as_dgC(x)
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = counts)
    )
  }

  cells <- colnames(sce)
  genes <- rownames(sce)

  merged_coldata <- .scprokar_align_data_frame(
    if (!is.null(cell_metadata)) cell_metadata else as.data.frame(SummarizedExperiment::colData(sce)),
    ids = cells,
    what = "cell_metadata"
  )
  merged_rowdata <- .scprokar_align_data_frame(
    if (!is.null(feature_metadata)) feature_metadata else as.data.frame(SummarizedExperiment::rowData(sce)),
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

  meta <- .scprokar_get_metadata(sce)
  meta$package <- "SCProkaR"
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

  .scprokar_set_metadata(sce, meta)
}

#' Merge multiple SCProkaR or SingleCellExperiment objects
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
MergeBacObjects <- function(..., objects = NULL, gene_mode = c("intersect", "union")) {
  gene_mode <- match.arg(gene_mode)
  dots <- list(...)
  if (!is.null(objects)) {
    dots <- c(dots, objects)
  }
  dots <- Filter(Negate(is.null), dots)

  if (length(dots) < 2) {
    stop("Provide at least two SingleCellExperiment objects to merge.", call. = FALSE)
  }
  if (!all(vapply(dots, methods::is, logical(1), "SingleCellExperiment"))) {
    stop("All inputs to MergeBacObjects() must be SingleCellExperiment objects.", call. = FALSE)
  }

  dots <- .scprokar_prepare_objects_for_merge(dots)

  gene_sets <- lapply(dots, rownames)
  genes <- if (gene_mode == "intersect") {
    Reduce(intersect, gene_sets)
  } else {
    Reduce(union, gene_sets)
  }
  if (length(genes) == 0) {
    stop("No genes remained after applying `gene_mode = \"", gene_mode, "\"`.", call. = FALSE)
  }

  counts_list <- lapply(dots, function(sce) {
    counts <- SummarizedExperiment::assay(sce, "counts")
    .scprokar_expand_counts(counts, genes)
  })
  merged_counts <- do.call(Matrix::cbind2, counts_list)
  merged_counts <- .scprokar_as_dgC(merged_counts)

  coldata_list <- lapply(dots, function(sce) {
    as.data.frame(SummarizedExperiment::colData(sce))
  })
  merged_coldata <- do.call(rbind, coldata_list)
  rownames(merged_coldata) <- unlist(lapply(dots, colnames), use.names = FALSE)

  rowdata_list <- lapply(dots, function(sce) {
    df <- as.data.frame(SummarizedExperiment::rowData(sce))
    df <- df[match(rownames(sce), rownames(df)), , drop = FALSE]
    rownames(df) <- rownames(sce)
    df
  })
  merged_rowdata <- .scprokar_merge_rowdata(rowdata_list, genes)

  merged <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = merged_counts),
    colData = S4Vectors::DataFrame(merged_coldata),
    rowData = S4Vectors::DataFrame(merged_rowdata)
  )

  if (all(vapply(dots, function(sce) "logcounts" %in% SummarizedExperiment::assayNames(sce), logical(1)))) {
    logcounts_list <- lapply(dots, function(sce) {
      logcounts <- SummarizedExperiment::assay(sce, "logcounts")
      .scprokar_expand_counts(logcounts, genes)
    })
    SummarizedExperiment::assay(merged, "logcounts") <- .scprokar_as_dgC(do.call(Matrix::cbind2, logcounts_list))
  }

  merged_reduction <- .scprokar_merge_reduced_dims2(dots)
  for (reduction_name in names(merged_reduction$kept)) {
    reduction_matrix <- merged_reduction$kept[[reduction_name]]
    SingleCellExperiment::reducedDim(merged, reduction_name) <- reduction_matrix
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
.scprokar_from_seurat <- function(x, seurat_assay = NULL, seurat_layer = "counts", transfer_reductions = TRUE) {
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
    colData = S4Vectors::DataFrame(.scprokar_align_data_frame(cell_metadata, colnames(counts), "cell_metadata")),
    rowData = S4Vectors::DataFrame(feature_metadata)
  )

  normalized <- tryCatch(
    .scprokar_extract_seurat_layer(x, assay = seurat_assay, layer = "data"),
    error = function(e) NULL
  )
  if (!is.null(normalized) &&
      nrow(normalized) == nrow(sce) &&
      ncol(normalized) == ncol(sce)) {
    SummarizedExperiment::assay(sce, "logcounts") <- .scprokar_as_dgC(normalized[rownames(sce), colnames(sce), drop = FALSE])
  }

  reduction_names <- character(0)
  if (isTRUE(transfer_reductions)) {
    reduction_names <- tryCatch(
      names(methods::slot(x, "reductions")),
      error = function(e) character(0)
    )
    for (reduction_name in reduction_names) {
      emb <- tryCatch(
        SeuratObject::Embeddings(object = x, reduction = reduction_name),
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
      SeuratObject::GetAssayData(object = x, assay = assay, layer = layer),
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
.scprokar_expand_counts <- function(mat, genes) {
  out <- Matrix::Matrix(0, nrow = length(genes), ncol = ncol(mat), sparse = TRUE)
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
  all_columns <- unique(unlist(lapply(rowdata_list, colnames), use.names = FALSE))
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
.scprokar_merge_reduced_dims <- function(objects) {
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
    ncell <- vapply(reductions, nrow, integer(1))
    if (length(unique(ndim)) != 1L) {
      skipped[[reduction_name]] <- sprintf(
        "dimension mismatch (ncol: %s)",
        paste(ndim, collapse = ",")
      )
      next
    }
    if (length(unique(ncell)) != 1L || any(ncell != vapply(objects, ncol, integer(1)))) {
      skipped[[reduction_name]] <- "cell-count mismatch between embedding rows and object columns"
      next
    }
    merged_matrix <- do.call(rbind, reductions)
    rownames(merged_matrix) <- cell_ids
    kept[[reduction_name]] <- merged_matrix
  }

  list(kept = kept, skipped = skipped, n_objects = length(objects))
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
      skipped[[reduction_name]] <- "unable to align reduced dimensions by cell IDs"
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
  for (candidate in c(cols$sample_col, cols$batch_col, "sample_id", "batch")) {
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
