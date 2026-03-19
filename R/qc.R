#' Compute bacterial QC metrics
#'
#' Adds Seurat-style per-cell feature fractions to a `SingleCellExperiment`,
#' including rRNA and ribosomal-protein content.
#'
#' @param sce A `SingleCellExperiment`.
#' @param rrna_pattern Regular expression used to identify rRNA features.
#' @param ribo_pattern Regular expression used to identify ribosomal protein
#'   features.
#' @param gene_class_col Optional `rowData` column describing feature classes.
#'   When present, values matching `"rrna"` are used for rRNA metrics and values
#'   matching `"ribo"` or `"ribosomal_protein"` are used for ribosomal protein
#'   metrics.
#' @param store If `TRUE`, store metrics in `colData(sce)` and return the
#'   modified object. If `FALSE`, return the QC table only.
#'
#' @return A `SingleCellExperiment` when `store = TRUE`, otherwise a data frame
#'   of QC metrics.
#' @export
RunBacQC <- function(
    sce,
    rrna_pattern = "^(rrs|rrl|rrf)",
    ribo_pattern = "^(rpl|rps)",
    gene_class_col = NULL,
    store = TRUE
) {
  counts <- SummarizedExperiment::assay(sce, "counts")
  if (is.null(counts)) {
    stop("No `counts` assay found. Run CreateBacObject() first.", call. = FALSE)
  }

  row_data <- as.data.frame(SummarizedExperiment::rowData(sce))
  genes <- rownames(sce)
  rrna_index <- grepl(rrna_pattern, genes, ignore.case = TRUE)
  ribo_index <- grepl(ribo_pattern, genes, ignore.case = TRUE)

  if (!is.null(gene_class_col)) {
    if (!gene_class_col %in% colnames(row_data)) {
      stop("`gene_class_col` was not found in rowData(sce).", call. = FALSE)
    }
    classes <- tolower(as.character(row_data[[gene_class_col]]))
    rrna_index <- classes %in% c("rrna", "r_rna", "ribosomal_rna")
    ribo_index <- classes %in% c("ribo", "ribosomal_protein", "ribo_protein")
  }

  total_counts <- Matrix::colSums(counts)
  detected_features <- Matrix::colSums(counts > 0)
  rrna_stats <- .scprokar_fraction_from_features(counts, rrna_index)
  ribo_stats <- .scprokar_fraction_from_features(counts, ribo_index)

  qc_df <- data.frame(
    total_counts = as.numeric(total_counts),
    detected_features = as.numeric(detected_features),
    rrna_counts = rrna_stats$selected,
    rrna_fraction = rrna_stats$fraction,
    pct_rrna = rrna_stats$fraction * 100,
    ribo_counts = ribo_stats$selected,
    ribo_fraction = ribo_stats$fraction,
    pct_ribo = ribo_stats$fraction * 100,
    row.names = colnames(sce)
  )

  if (!isTRUE(store)) {
    return(qc_df)
  }

  col_data <- as.data.frame(SummarizedExperiment::colData(sce))
  for (nm in colnames(qc_df)) {
    col_data[[nm]] <- qc_df[[nm]]
  }
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(col_data)

  .scprokar_store_step(
    sce,
    "qc",
    list(
      rrna_pattern = rrna_pattern,
      ribo_pattern = ribo_pattern,
      gene_class_col = gene_class_col,
      metrics = colnames(qc_df)
    )
  )
}

#' Filter cells using bacterial QC thresholds
#'
#' Applies simple threshold-based filtering to QC metrics produced by
#' `RunBacQC()`.
#'
#' @param sce A `SingleCellExperiment` with QC metrics already present.
#' @param min_counts Minimum allowed total counts.
#' @param min_features Minimum allowed number of detected features.
#' @param max_rrna_fraction Maximum allowed rRNA fraction.
#' @param custom_filter Optional custom filter. May be a logical vector, a
#'   column name in `colData(sce)`, or a function returning a logical vector.
#'
#' @return A filtered `SingleCellExperiment`.
#' @export
FilterBacCells <- function(
    sce,
    min_counts = NULL,
    min_features = NULL,
    max_rrna_fraction = NULL,
    custom_filter = NULL
) {
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  required <- c("total_counts", "detected_features", "rrna_fraction")
  missing_cols <- setdiff(required, colnames(cd))
  if (length(missing_cols) > 0) {
    stop("Run RunBacQC() before filtering cells.", call. = FALSE)
  }

  keep <- rep(TRUE, nrow(cd))
  if (!is.null(min_counts)) {
    keep <- keep & (cd$total_counts >= min_counts)
  }
  if (!is.null(min_features)) {
    keep <- keep & (cd$detected_features >= min_features)
  }
  if (!is.null(max_rrna_fraction)) {
    keep <- keep & (cd$rrna_fraction <= max_rrna_fraction)
  }

  if (!is.null(custom_filter)) {
    if (is.function(custom_filter)) {
      custom_keep <- custom_filter(cd)
    } else if (is.character(custom_filter) && length(custom_filter) == 1L) {
      custom_keep <- cd[[custom_filter]]
    } else {
      custom_keep <- custom_filter
    }
    if (!is.logical(custom_keep) || length(custom_keep) != nrow(cd)) {
      stop("`custom_filter` must resolve to one logical value per cell.", call. = FALSE)
    }
    keep <- keep & custom_keep
  }

  filtered <- sce[, keep]
  .scprokar_store_step(
    filtered,
    "filtering",
    list(
      min_counts = min_counts,
      min_features = min_features,
      max_rrna_fraction = max_rrna_fraction,
      kept_cells = sum(keep),
      removed_cells = sum(!keep)
    )
  )
}
