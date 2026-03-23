#' Aggregate single-cell data into pseudobulk samples
#'
#' Summarize single-cell counts into pseudobulk samples defined by one or more
#' metadata columns. This is a common first step before sample-level
#' differential expression analysis with tools such as edgeR.
#'
#' @param sce A `SingleCellExperiment`.
#' @param sample_cols Character vector of `colData(sce)` columns that together
#'   define each pseudobulk sample, for example `c("sample_id", "timepoint",
#'   "treatment", "cluster")`.
#' @param assay_name Assay containing the input counts.
#' @param aggregation Aggregation method. Choices are:
#' - `"sum"`: sum counts across cells in each pseudobulk sample. This is the
#'   recommended default for count-based edgeR workflows.
#' - `"mean"`: average expression across cells in each pseudobulk sample. This
#'   can be useful for exploratory summaries but is not ideal for count-based
#'   NB models.
#' @param min_cells Minimum number of cells required for a pseudobulk sample to
#'   be retained.
#' @param add_logcounts Logical indicating whether to add a `logcounts` assay
#'   using `edgeR::cpm(..., log = TRUE)`.
#' @param prior_count Prior count passed to `edgeR::cpm()` when `add_logcounts =
#'   TRUE`.
#' @param sample_prefix Prefix for the generated pseudobulk column names.
#'
#' @return A pseudobulk `SingleCellExperiment` with one column per aggregated
#'   sample. The output `colData` includes the grouping metadata, `ncells`, and
#'   `lib.size`.
#' @export
aggregate_pseudobulk <- function(
    sce,
    sample_cols,
    assay_name = "counts",
    aggregation = c("sum", "mean"),
    min_cells = 1L,
    add_logcounts = TRUE,
    prior_count = 2,
    sample_prefix = "pb_") {
  aggregation <- match.arg(aggregation)

  if (!methods::is(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment.", call. = FALSE)
  }

  if (!assay_name %in% SummarizedExperiment::assayNames(sce)) {
    stop("Assay `", assay_name, "` is not present in `sce`.", call. = FALSE)
  }

  if (length(sample_cols) == 0L || !all(sample_cols %in% colnames(SummarizedExperiment::colData(sce)))) {
    stop("`sample_cols` must all be present in `colData(sce)`.", call. = FALSE)
  }

  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  keep_cells <- stats::complete.cases(meta[, sample_cols, drop = FALSE])
  if (!any(keep_cells)) {
    stop("No cells have complete values across `sample_cols`.", call. = FALSE)
  }

  meta <- meta[keep_cells, , drop = FALSE]
  counts <- SummarizedExperiment::assay(sce, assay_name)[, keep_cells, drop = FALSE]

  group_df <- meta[, sample_cols, drop = FALSE]
  group_factor <- interaction(group_df, drop = TRUE, lex.order = TRUE, sep = "||")
  group_levels <- levels(group_factor)
  ncells <- as.integer(table(group_factor))
  keep_groups <- ncells >= as.integer(min_cells)

  if (!any(keep_groups)) {
    stop("No pseudobulk samples remain after applying `min_cells`.", call. = FALSE)
  }

  membership <- Matrix::sparseMatrix(
    i = seq_along(group_factor),
    j = as.integer(group_factor),
    x = 1,
    dims = c(length(group_factor), length(group_levels))
  )

  agg_counts <- counts %*% membership
  if (aggregation == "mean") {
    agg_counts <- t(t(agg_counts) / pmax(ncells, 1))
  }

  agg_counts <- agg_counts[, keep_groups, drop = FALSE]
  ncells <- ncells[keep_groups]
  kept_levels <- group_levels[keep_groups]

  sample_meta <- unique(data.frame(.group = as.character(group_factor), group_df, stringsAsFactors = FALSE))
  sample_meta <- sample_meta[match(kept_levels, sample_meta$.group), , drop = FALSE]
  rownames(sample_meta) <- paste0(sample_prefix, seq_len(nrow(sample_meta)))
  sample_meta$.group <- NULL
  sample_meta$ncells <- ncells
  sample_meta$lib.size <- Matrix::colSums(agg_counts)

  colnames(agg_counts) <- rownames(sample_meta)
  rownames(agg_counts) <- rownames(sce)

  pb <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = agg_counts)
  )
  SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(sample_meta)

  if (isTRUE(add_logcounts)) {
    SummarizedExperiment::assay(pb, "logcounts") <- edgeR::cpm(
      SummarizedExperiment::assay(pb, "counts"),
      log = TRUE,
      prior.count = prior_count
    )
  }

  S4Vectors::metadata(pb)$pseudobulk <- list(
    sample_cols = sample_cols,
    aggregation = aggregation,
    assay_name = assay_name,
    min_cells = min_cells
  )

  pb
}


#' Filter pseudobulk samples by size
#'
#' Filter a pseudobulk `SingleCellExperiment` using minimum numbers of cells and
#' library size.
#'
#' @param pb A pseudobulk `SingleCellExperiment`, typically returned by
#'   `aggregate_pseudobulk()`.
#' @param min_cells Minimum number of cells required.
#' @param min_lib_size Minimum library size required.
#'
#' @return A filtered pseudobulk `SingleCellExperiment`.
#' @export
filter_pseudobulk_samples <- function(
    pb,
    min_cells = 10L,
    min_lib_size = NULL) {
  if (!methods::is(pb, "SingleCellExperiment")) {
    stop("`pb` must be a SingleCellExperiment.", call. = FALSE)
  }

  meta <- as.data.frame(SummarizedExperiment::colData(pb))
  keep <- rep(TRUE, ncol(pb))

  if ("ncells" %in% colnames(meta) && !is.null(min_cells)) {
    keep <- keep & meta$ncells >= as.integer(min_cells)
  }

  if (is.null(min_lib_size) && "lib.size" %in% colnames(meta)) {
    min_lib_size <- 0
  }
  if (!is.null(min_lib_size)) {
    lib_size <- if ("lib.size" %in% colnames(meta)) {
      meta$lib.size
    } else {
      Matrix::colSums(SummarizedExperiment::assay(pb, "counts"))
    }
    keep <- keep & lib_size >= min_lib_size
  }

  pb[, keep, drop = FALSE]
}


#' Normalize pseudobulk samples for exploratory analysis
#'
#' Compute edgeR normalization factors and add a `logcounts` assay for plotting
#' and exploratory analysis.
#'
#' @param pb A pseudobulk `SingleCellExperiment`.
#' @param assay_name Assay containing pseudobulk counts.
#' @param method Normalization method. Choices are:
#' - `"TMM"`: trimmed mean of M values normalization via
#'   `edgeR::calcNormFactors(..., method = "TMM")`.
#' - `"upperquartile"`: upper-quartile normalization via edgeR.
#' - `"none"`: do not apply normalization factors; still compute `logcounts`
#'   from the raw library sizes.
#' @param output_assay Name of the output log-scale assay.
#' @param prior_count Prior count passed to `edgeR::cpm(..., log = TRUE)`.
#'
#' @return The input pseudobulk object with updated `colData` normalization
#'   fields and a log-scale assay.
#' @export
normalize_pseudobulk <- function(
    pb,
    assay_name = "counts",
    method = c("TMM", "upperquartile", "none"),
    output_assay = "logcounts",
    prior_count = 2) {
  method <- match.arg(method)

  if (!methods::is(pb, "SingleCellExperiment")) {
    stop("`pb` must be a SingleCellExperiment.", call. = FALSE)
  }

  if (!assay_name %in% SummarizedExperiment::assayNames(pb)) {
    stop("Assay `", assay_name, "` is not present in `pb`.", call. = FALSE)
  }

  dge <- edgeR::DGEList(counts = SummarizedExperiment::assay(pb, assay_name))
  if (method != "none") {
    dge <- edgeR::calcNormFactors(dge, method = method)
  } else {
    dge$samples$norm.factors <- rep(1, ncol(pb))
  }

  SummarizedExperiment::colData(pb)$lib.size <- dge$samples$lib.size
  SummarizedExperiment::colData(pb)$norm.factors <- dge$samples$norm.factors
  SummarizedExperiment::assay(pb, output_assay) <- edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = prior_count
  )

  pb
}
