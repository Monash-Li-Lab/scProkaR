#' Aggregate cells into pseudobulk profiles
#'
#' Sums or averages counts across cell groups to prepare a pseudobulk analysis
#' object for downstream edgeR modeling.
#'
#' @param sce A `SingleCellExperiment`.
#' @param sample_col Column containing sample identifiers.
#' @param group_cols Additional grouping columns used to define pseudobulk
#'   profiles.
#' @param assay_name Assay to aggregate.
#' @param fun Aggregation function, either `"sum"` or `"mean"`.
#'
#' @return A pseudobulk object of class `SCProkaRPseudobulk`.
#' @export
AggregatePseudobulk <- function(
    sce,
    sample_col,
    group_cols,
    assay_name = "counts",
    fun = "sum"
) {
  fun <- match.arg(fun, c("sum", "mean"))
  .scprokar_match_columns(sce, c(sample_col, group_cols), label = "grouping")

  counts <- SummarizedExperiment::assay(sce, assay_name)
  sample_df <- as.data.frame(SummarizedExperiment::colData(sce)[, c(sample_col, group_cols), drop = FALSE])
  group_id <- do.call(interaction, c(sample_df, list(drop = TRUE, lex.order = TRUE)))
  levels_id <- levels(group_id)

  aggregated <- lapply(levels_id, function(id) {
    idx <- which(group_id == id)
    if (fun == "sum") {
      Matrix::rowSums(counts[, idx, drop = FALSE])
    } else {
      Matrix::rowMeans(counts[, idx, drop = FALSE])
    }
  })
  agg_matrix <- do.call(cbind, aggregated)
  colnames(agg_matrix) <- levels_id
  rownames(agg_matrix) <- rownames(sce)

  samples <- unique(data.frame(group_id = as.character(group_id), sample_df, stringsAsFactors = FALSE))
  rownames(samples) <- samples$group_id
  samples <- samples[colnames(agg_matrix), , drop = FALSE]
  samples$n_cells <- as.integer(table(group_id)[colnames(agg_matrix)])

  out <- list(
    counts = .scprokar_as_dgC(agg_matrix),
    samples = samples,
    sample_col = sample_col,
    group_cols = group_cols,
    assay_name = assay_name,
    aggregation = fun
  )
  class(out) <- "SCProkaRPseudobulk"
  out
}

#' Run pseudobulk differential expression with edgeR
#'
#' Fits an edgeR quasi-likelihood model on pseudobulk counts and returns a
#' standardized result table.
#'
#' @param pb A pseudobulk object from `AggregatePseudobulk()`.
#' @param design A formula or design matrix.
#' @param contrast A coefficient name, contrast vector, or one-column contrast
#'   matrix.
#' @param method edgeR backend name. Both `"edgeR"` and `"QLF"` use the
#'   quasi-likelihood workflow.
#'
#' @return A list containing the fitted model and a result table.
#' @export
RunPseudobulkDE <- function(pb, design, contrast, method = c("edgeR", "QLF")) {
  method <- match.arg(method)
  .scprokar_require("edgeR", "pseudobulk differential expression")

  if (!inherits(pb, "SCProkaRPseudobulk")) {
    stop("`pb` must be created by AggregatePseudobulk().", call. = FALSE)
  }

  design_matrix <- if (inherits(design, "formula")) {
    stats::model.matrix(design, data = pb$samples)
  } else {
    as.matrix(design)
  }

  if (is.character(contrast) && length(contrast) == 1L) {
    if (!contrast %in% colnames(design_matrix)) {
      stop(
        "`contrast` must match a column in the design matrix. Available columns: ",
        paste(colnames(design_matrix), collapse = ", "),
        call. = FALSE
      )
    }
  } else if (is.numeric(contrast)) {
    if (length(contrast) != ncol(design_matrix)) {
      stop("Numeric `contrast` must have one value per design-matrix column.", call. = FALSE)
    }
  } else if (is.matrix(contrast)) {
    if (ncol(contrast) != 1L || nrow(contrast) != ncol(design_matrix)) {
      stop("Matrix `contrast` must be one column with one row per design-matrix column.", call. = FALSE)
    }
  } else {
    stop("`contrast` must be a coefficient name, numeric vector, or one-column matrix.", call. = FALSE)
  }

  y <- edgeR::DGEList(counts = pb$counts, samples = pb$samples)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design = design_matrix)
  fit <- edgeR::glmQLFit(y, design = design_matrix)

  test <- if (is.character(contrast) && length(contrast) == 1L) {
    edgeR::glmQLFTest(fit, coef = contrast)
  } else if (is.numeric(contrast)) {
    edgeR::glmQLFTest(fit, contrast = contrast)
  } else if (is.matrix(contrast)) {
    edgeR::glmQLFTest(fit, contrast = contrast[, 1])
  }

  results <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  results$gene_id <- rownames(results)
  results <- results[, c("gene_id", setdiff(colnames(results), "gene_id")), drop = FALSE]

  list(
    method = method,
    fit = fit,
    test = test,
    design_matrix = design_matrix,
    contrast = contrast,
    results = results
  )
}

#' Run time-series or lineage-aware differential expression
#'
#' Fits a `tradeSeq` generalized additive model on raw counts using either
#' pseudotime inferred by `RunBacPAGA()` or a user-provided `time_col`.
#'
#' @param sce A `SingleCellExperiment`.
#' @param cluster_col Optional cluster column for downstream reporting.
#' @param time_col Column providing numeric ordering when trajectory output is
#'   unavailable.
#' @param condition_col Optional condition column for condition-aware tests.
#' @param lineage_col Optional lineage column used to create lineage weights.
#' @param method Time-series backend. Currently only `"tradeSeq"` is supported.
#'
#' @return A list with the fitted GAM and a standardized result table.
#' @export
RunTimeSeriesDE <- function(
    sce,
    cluster_col = NULL,
    time_col,
    condition_col = NULL,
    lineage_col = NULL,
    method = "tradeSeq"
) {
  if (!identical(method, "tradeSeq")) {
    stop("Only `tradeSeq` is supported in v1.", call. = FALSE)
  }
  .scprokar_require("tradeSeq", "time-series differential expression")

  col_data <- as.data.frame(SummarizedExperiment::colData(sce))
  if (!time_col %in% colnames(col_data) && !"pseudotime" %in% colnames(col_data)) {
    stop("`time_col` was not found and no pseudotime is stored in `colData(sce)`.", call. = FALSE)
  }

  if (!is.null(cluster_col)) {
    .scprokar_match_columns(sce, cluster_col, label = "cluster")
  }
  if (!is.null(condition_col)) {
    .scprokar_match_columns(sce, condition_col, label = "condition")
  }
  if (!is.null(lineage_col)) {
    .scprokar_match_columns(sce, lineage_col, label = "lineage")
  }

  counts <- SummarizedExperiment::assay(sce, "counts")
  pseudotime_vector <- if ("pseudotime" %in% colnames(col_data)) {
    as.numeric(col_data$pseudotime)
  } else {
    as.numeric(col_data[[time_col]])
  }

  if (!is.null(lineage_col)) {
    lineage <- factor(col_data[[lineage_col]])
    pseudotime <- matrix(rep(pseudotime_vector, length(levels(lineage))), ncol = length(levels(lineage)))
    colnames(pseudotime) <- levels(lineage)
    cell_weights <- stats::model.matrix(~ lineage - 1)
    colnames(cell_weights) <- levels(lineage)
  } else {
    pseudotime <- matrix(pseudotime_vector, ncol = 1)
    colnames(pseudotime) <- "Lineage1"
    cell_weights <- matrix(1, nrow = length(pseudotime_vector), ncol = 1)
    colnames(cell_weights) <- "Lineage1"
  }

  fit_args <- list(
    counts = counts,
    pseudotime = pseudotime,
    cellWeights = cell_weights,
    verbose = FALSE
  )
  if (!is.null(condition_col)) {
    fit_args$conditions <- factor(col_data[[condition_col]])
  }
  gam <- do.call(tradeSeq::fitGAM, fit_args)

  result_table <- if (!is.null(condition_col)) {
    tradeSeq::conditionTest(gam)
  } else {
    tradeSeq::associationTest(gam)
  }

  result_table$gene_id <- rownames(result_table)
  result_table <- result_table[, c("gene_id", setdiff(colnames(result_table), "gene_id")), drop = FALSE]

  list(
    method = method,
    gam = gam,
    cluster_col = cluster_col,
    condition_col = condition_col,
    lineage_col = lineage_col,
    results = result_table
  )
}
