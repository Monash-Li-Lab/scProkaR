# Internal helpers for the edgeR workflows.

.most_frequent_value <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(NA)
  }
  tab <- table(x)
  names(tab)[which.max(tab)]
}


.sanitize_levels <- function(x) {
  stats::setNames(make.names(x), x)
}


.make_default_reference <- function(values) {
  if (is.numeric(values) || is.integer(values)) {
    unique_vals <- sort(unique(as.numeric(values)))
    return(as.character(unique_vals[1]))
  }

  if (is.factor(values)) {
    return(as.character(levels(values)[1]))
  }

  unique_vals <- unique(as.character(values))
  suppressWarnings(numeric_vals <- as.numeric(unique_vals))
  if (!anyNA(numeric_vals)) {
    return(as.character(unique_vals[order(numeric_vals)][1]))
  }
  sort(unique_vals)[1]
}


.build_group_factor <- function(meta, group_col, combine_group_cols = NULL) {
  if (is.null(combine_group_cols)) {
    values <- meta[[group_col]]
    if (is.numeric(values) || is.integer(values)) {
      ordered_levels <- as.character(sort(unique(values)))
      return(factor(as.character(values), levels = ordered_levels))
    }
    if (is.factor(values)) {
      return(factor(as.character(values), levels = levels(values)))
    }
    values_chr <- as.character(values)
    suppressWarnings(values_num <- as.numeric(values_chr))
    if (!anyNA(values_num)) {
      ordered_levels <- unique(values_chr[order(values_num)])
      return(factor(values_chr, levels = ordered_levels))
    }
    return(factor(values_chr, levels = sort(unique(values_chr))))
  }

  combined <- interaction(meta[, combine_group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE, sep = "__")
  factor(as.character(combined), levels = levels(combined))
}


.make_pairwise_pairs <- function(group_levels, contrast_type, reference_level = NULL, manual_pairs = NULL) {
  contrast_type <- match.arg(contrast_type, c("reference", "all", "manual"))

  if (contrast_type == "manual") {
    if (is.null(manual_pairs)) {
      stop("`manual_pairs` must be supplied when `contrast_type = \"manual\"`.", call. = FALSE)
    }
    if (!is.matrix(manual_pairs) && !is.data.frame(manual_pairs)) {
      stop("`manual_pairs` must be a two-column matrix or data.frame.", call. = FALSE)
    }
    manual_pairs <- as.data.frame(manual_pairs, stringsAsFactors = FALSE)
    if (ncol(manual_pairs) != 2L) {
      stop("`manual_pairs` must have exactly two columns.", call. = FALSE)
    }
    colnames(manual_pairs) <- c("group1", "group2")
    return(manual_pairs)
  }

  if (contrast_type == "reference") {
    if (is.null(reference_level)) {
      reference_level <- group_levels[1]
    }
    if (!reference_level %in% group_levels) {
      stop("`reference_level` is not among the group levels.", call. = FALSE)
    }
    out <- data.frame(
      group1 = setdiff(group_levels, reference_level),
      group2 = reference_level,
      stringsAsFactors = FALSE
    )
    return(out)
  }

  combn_out <- utils::combn(group_levels, 2)
  data.frame(
    group1 = combn_out[2, ],
    group2 = combn_out[1, ],
    stringsAsFactors = FALSE
  )
}


.prepare_edger_dge <- function(counts, meta, normalization = "TMM") {
  dge <- edgeR::DGEList(counts = as.matrix(counts), samples = meta)
  if (!is.null(normalization) && normalization != "none") {
    dge <- edgeR::calcNormFactors(dge, method = normalization)
  } else {
    dge$samples$norm.factors <- rep(1, ncol(counts))
  }
  dge
}


.resolve_pairwise_de_table <- function(de_result, contrast) {
  if (is.null(de_result$tables) || !length(de_result$tables)) {
    stop("`de_result` does not contain pairwise DE tables.", call. = FALSE)
  }

  if (is.null(contrast)) {
    contrast <- names(de_result$tables)[1]
  }

  if (!contrast %in% names(de_result$tables)) {
    stop("Contrast `", contrast, "` is not present in `de_result$tables`.", call. = FALSE)
  }

  tab <- de_result$tables[[contrast]]
  if (!all(c("gene", "logFC", "FDR") %in% colnames(tab))) {
    stop("The selected DE table does not contain the required columns.", call. = FALSE)
  }

  tab
}


#' Run pairwise pseudobulk differential expression with edgeR
#'
#' Fit a quasi-likelihood negative-binomial model with edgeR and test pairwise
#' contrasts between groups such as timepoints, treatments, or composite
#' treatment-timepoint groups.
#'
#' @param pb A pseudobulk `SingleCellExperiment`.
#' @param assay_name Assay containing pseudobulk counts.
#' @param group_col Sample metadata column defining the primary grouping, for
#'   example `"timepoint"`.
#' @param covariates Optional metadata columns to include as additive covariates
#'   in the design matrix, for example treatment, batch, or donor.
#' @param combine_group_cols Optional metadata columns to combine into a single
#'   composite group via `interaction()`, for example
#'   `c("treatment", "timepoint")`.
#' @param contrast_type Contrast construction strategy. Choices are:
#' - `"reference"`: compare every non-reference group against one reference
#'   level. This is the default.
#' - `"all"`: test all pairwise group contrasts.
#' - `"manual"`: use `manual_pairs` exactly as supplied.
#' @param reference_level Optional reference level used when
#'   `contrast_type = "reference"`. By default the earliest or first group is
#'   used.
#' @param manual_pairs Optional two-column matrix or data frame of manual
#'   contrasts with `group1` and `group2` values. Each test is `group1 - group2`.
#' @param normalization Library-size normalization method. Choices are:
#' - `"TMM"`: edgeR TMM normalization.
#' - `"upperquartile"`: edgeR upper-quartile normalization.
#' - `"none"`: no normalization factors beyond raw library sizes.
#' @param filter_by_expr Logical indicating whether to run
#'   `edgeR::filterByExpr()`.
#' @param robust Logical passed to `edgeR::glmQLFit()`.
#' @param abundance_trend Logical passed to `edgeR::glmQLFit()`.
#' @param fdr_cutoff FDR cutoff used to mark significant genes.
#'
#' @return A list containing the fitted edgeR objects, contrast matrix, result
#'   tables, and analysis settings.
#' @export
run_edger_pairwise_de <- function(
    pb,
    assay_name = "counts",
    group_col = "timepoint",
    covariates = NULL,
    combine_group_cols = NULL,
    contrast_type = c("reference", "all", "manual"),
    reference_level = NULL,
    manual_pairs = NULL,
    normalization = c("TMM", "upperquartile", "none"),
    filter_by_expr = TRUE,
    robust = TRUE,
    abundance_trend = FALSE,
    fdr_cutoff = 0.05) {
  contrast_type <- match.arg(contrast_type)
  normalization <- match.arg(normalization)

  if (!methods::is(pb, "SingleCellExperiment")) {
    stop("`pb` must be a SingleCellExperiment.", call. = FALSE)
  }
  if (!assay_name %in% SummarizedExperiment::assayNames(pb)) {
    stop("Assay `", assay_name, "` is not present in `pb`.", call. = FALSE)
  }

  meta <- as.data.frame(SummarizedExperiment::colData(pb))
  needed_cols <- unique(c(group_col, combine_group_cols, covariates))
  if (!all(needed_cols %in% colnames(meta))) {
    stop("All grouping and covariate columns must exist in `colData(pb)`.", call. = FALSE)
  }

  keep <- stats::complete.cases(meta[, needed_cols, drop = FALSE])
  if (!any(keep)) {
    stop("No pseudobulk samples have complete values for the requested model.", call. = FALSE)
  }

  meta <- meta[keep, , drop = FALSE]
  counts <- SummarizedExperiment::assay(pb, assay_name)[, keep, drop = FALSE]

  group_factor <- .build_group_factor(meta, group_col = group_col, combine_group_cols = combine_group_cols)
  group_levels <- levels(group_factor)
  if (length(group_levels) < 2L) {
    stop("At least two groups are required for pairwise DE.", call. = FALSE)
  }

  if (is.null(reference_level) && contrast_type == "reference") {
    reference_level <- .make_default_reference(meta[[if (is.null(combine_group_cols)) group_col else combine_group_cols[1]]])
    if (!reference_level %in% group_levels) {
      reference_level <- group_levels[1]
    }
  }

  design_df <- data.frame(group_factor = group_factor, stringsAsFactors = FALSE)
  if (!is.null(covariates) && length(covariates) > 0L) {
    design_df <- cbind(design_df, meta[, covariates, drop = FALSE])
  }

  design <- stats::model.matrix(~ 0 + ., data = design_df)
  safe_group_levels <- .sanitize_levels(group_levels)
  group_cols <- seq_along(group_levels)
  colnames(design)[group_cols] <- paste0("group__", unname(safe_group_levels))
  colnames(design) <- make.names(colnames(design))

  dge <- .prepare_edger_dge(counts = counts, meta = meta, normalization = normalization)
  if (isTRUE(filter_by_expr)) {
    keep_genes <- edgeR::filterByExpr(dge, design = design, group = group_factor)
    dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  }

  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = robust, abundance.trend = abundance_trend)

  pair_df <- .make_pairwise_pairs(
    group_levels = group_levels,
    contrast_type = contrast_type,
    reference_level = reference_level,
    manual_pairs = manual_pairs
  )

  contrast_names <- paste0(pair_df$group1, "_vs_", pair_df$group2)
  contrast_matrix <- matrix(
    0,
    nrow = ncol(design),
    ncol = nrow(pair_df),
    dimnames = list(colnames(design), contrast_names)
  )
  test_list <- vector("list", nrow(pair_df))
  table_list <- vector("list", nrow(pair_df))

  for (i in seq_len(nrow(pair_df))) {
    g1 <- pair_df$group1[i]
    g2 <- pair_df$group2[i]
    if (!g1 %in% group_levels || !g2 %in% group_levels) {
      stop("All manual contrast groups must appear in the design.", call. = FALSE)
    }

    contrast_name <- contrast_names[i]
    contrast_vec <- rep(0, ncol(design))
    contrast_vec[match(paste0("group__", safe_group_levels[[g1]]), colnames(design))] <- 1
    contrast_vec[match(paste0("group__", safe_group_levels[[g2]]), colnames(design))] <- -1
    contrast_matrix[, i] <- contrast_vec

    qlf <- edgeR::glmQLFTest(fit, contrast = contrast_vec)
    tab <- edgeR::topTags(qlf, n = Inf)$table
    tab$gene <- rownames(tab)
    tab$contrast <- contrast_name
    tab$group1 <- g1
    tab$group2 <- g2
    tab$significant <- tab$FDR <= fdr_cutoff

    test_list[[i]] <- qlf
    table_list[[i]] <- tab
  }

  names(test_list) <- colnames(contrast_matrix)
  names(table_list) <- colnames(contrast_matrix)

  list(
    dge = dge,
    design = design,
    fit = fit,
    group_factor = group_factor,
    contrast_matrix = contrast_matrix,
    tests = test_list,
    tables = table_list,
    combined_table = do.call(rbind, table_list),
    settings = list(
      assay_name = assay_name,
      group_col = group_col,
      covariates = covariates,
      combine_group_cols = combine_group_cols,
      contrast_type = contrast_type,
      reference_level = reference_level,
      normalization = normalization,
      filter_by_expr = filter_by_expr,
      robust = robust,
      abundance_trend = abundance_trend,
      fdr_cutoff = fdr_cutoff
    )
  )
}


#' Plot a volcano plot for a pairwise edgeR contrast
#'
#' Visualize one pairwise contrast from [run_edger_pairwise_de()] as a volcano
#' plot using log-fold change and FDR.
#'
#' @param de_result Result list returned by `run_edger_pairwise_de()`.
#' @param contrast Name of the contrast to plot. If `NULL`, the first available
#'   contrast is used.
#' @param fdr_cutoff FDR threshold used to define significant genes.
#' @param lfc_cutoff Absolute log-fold-change threshold used to highlight larger
#'   effects.
#' @param top_n_labels Number of top genes to label by FDR. Set to `0` to omit
#'   labels.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#'
#' @return A `ggplot2` object.
#' @export
plot_pairwise_de_volcano <- function(
    de_result,
    contrast = NULL,
    fdr_cutoff = 0.05,
    lfc_cutoff = 0,
    top_n_labels = 10L,
    point_size = 1.5,
    point_alpha = 0.75) {
  tab <- .resolve_pairwise_de_table(de_result, contrast = contrast)
  contrast_name <- unique(tab$contrast)[1]

  tab$neg_log10_fdr <- -log10(pmax(tab$FDR, .Machine$double.xmin))
  tab$de_class <- ifelse(
    tab$FDR <= fdr_cutoff & tab$logFC >= lfc_cutoff,
    "Up",
    ifelse(
      tab$FDR <= fdr_cutoff & tab$logFC <= -lfc_cutoff,
      "Down",
      "NotSig"
    )
  )
  tab$de_class <- factor(tab$de_class, levels = c("Down", "NotSig", "Up"))

  p <- ggplot2::ggplot(tab, ggplot2::aes(x = logFC, y = neg_log10_fdr, colour = de_class)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2, colour = "grey60") +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, colour = "grey60") +
    ggplot2::scale_colour_manual(values = c(Down = "#2c7bb6", NotSig = "grey70", Up = "#d7191c")) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste("Volcano plot:", contrast_name),
      x = "logFC",
      y = expression(-log[10](FDR)),
      colour = "status"
    )

  top_n_labels <- as.integer(top_n_labels)
  if (!is.na(top_n_labels) && top_n_labels > 0L) {
    label_df <- tab[order(tab$FDR, -abs(tab$logFC)), c("gene", "logFC", "neg_log10_fdr"), drop = FALSE]
    label_df <- utils::head(label_df, top_n_labels)
    p <- p + ggplot2::geom_text(
      data = label_df,
      mapping = ggplot2::aes(x = logFC, y = neg_log10_fdr, label = gene),
      inherit.aes = FALSE,
      size = 3,
      vjust = -0.4,
      check_overlap = TRUE
    )
  }

  p
}


#' Plot an MA plot for a pairwise edgeR contrast
#'
#' Visualize one pairwise contrast from [run_edger_pairwise_de()] as an MA plot
#' using average abundance and log-fold change.
#'
#' @param de_result Result list returned by `run_edger_pairwise_de()`.
#' @param contrast Name of the contrast to plot. If `NULL`, the first available
#'   contrast is used.
#' @param fdr_cutoff FDR threshold used to define significant genes.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#'
#' @return A `ggplot2` object.
#' @export
plot_pairwise_de_ma <- function(
    de_result,
    contrast = NULL,
    fdr_cutoff = 0.05,
    point_size = 1.5,
    point_alpha = 0.75) {
  tab <- .resolve_pairwise_de_table(de_result, contrast = contrast)
  contrast_name <- unique(tab$contrast)[1]

  if (!"logCPM" %in% colnames(tab)) {
    stop("The selected DE table does not contain a `logCPM` column for MA plotting.", call. = FALSE)
  }

  tab$de_class <- ifelse(
    tab$FDR <= fdr_cutoff & tab$logFC > 0,
    "Up",
    ifelse(tab$FDR <= fdr_cutoff & tab$logFC < 0, "Down", "NotSig")
  )
  tab$de_class <- factor(tab$de_class, levels = c("Down", "NotSig", "Up"))

  ggplot2::ggplot(tab, ggplot2::aes(x = logCPM, y = logFC, colour = de_class)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::scale_colour_manual(values = c(Down = "#2c7bb6", NotSig = "grey70", Up = "#d7191c")) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste("MA plot:", contrast_name),
      x = "average logCPM",
      y = "logFC",
      colour = "status"
    )
}


.fit_one_spline_de <- function(
    counts,
    meta,
    time_col,
    covariates,
    normalization,
    df,
    filter_by_expr,
    robust,
    abundance_trend,
    fdr_cutoff,
    return_curve_for,
    curve_grid_length,
    label = "all") {
  time_numeric <- .coerce_time_to_numeric(meta[[time_col]])
  if (length(unique(time_numeric)) < 3L) {
    stop("Spline DE requires at least three distinct time values.", call. = FALSE)
  }

  effective_df <- min(as.integer(df), length(unique(time_numeric)) - 1L)
  if (effective_df < 1L) {
    stop("Spline DE requires at least one estimable spline degree of freedom.", call. = FALSE)
  }

  spline_basis <- splines::ns(time_numeric, df = effective_df)
  spline_df <- as.data.frame(spline_basis)
  colnames(spline_df) <- paste0("spline_", seq_len(ncol(spline_df)))

  design_df <- spline_df
  if (!is.null(covariates) && length(covariates) > 0L) {
    design_df <- cbind(design_df, meta[, covariates, drop = FALSE])
  }
  design_formula <- stats::as.formula(paste("~", paste(colnames(design_df), collapse = " + ")))
  design <- stats::model.matrix(design_formula, data = design_df)
  design <- as.matrix(design)

  dge <- .prepare_edger_dge(counts = counts, meta = meta, normalization = normalization)
  if (isTRUE(filter_by_expr)) {
    keep_genes <- edgeR::filterByExpr(dge, design = design)
    dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  }

  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = robust, abundance.trend = abundance_trend)
  spline_coef <- seq_len(ncol(spline_df)) + 1L
  qlf <- edgeR::glmQLFTest(fit, coef = spline_coef)

  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene <- rownames(tab)
  tab$condition <- label
  tab$significant <- tab$FDR <= fdr_cutoff

  if (return_curve_for == "significant") {
    curve_genes <- tab$gene[tab$significant]
  } else {
    curve_genes <- tab$gene
  }

  curve_df <- data.frame()
  if (length(curve_genes) > 0L) {
    time_grid <- seq(min(time_numeric), max(time_numeric), length.out = curve_grid_length)
    new_basis <- as.data.frame(stats::predict(spline_basis, newx = time_grid))
    colnames(new_basis) <- colnames(spline_df)

    if (!is.null(covariates) && length(covariates) > 0L) {
      template <- meta[rep(1, curve_grid_length), covariates, drop = FALSE]
      for (cov_name in covariates) {
        if (is.numeric(meta[[cov_name]]) || is.integer(meta[[cov_name]])) {
          template[[cov_name]] <- mean(meta[[cov_name]], na.rm = TRUE)
        } else {
          mode_value <- .most_frequent_value(meta[[cov_name]])
          if (is.factor(meta[[cov_name]])) {
            template[[cov_name]] <- factor(rep(mode_value, curve_grid_length), levels = levels(meta[[cov_name]]))
          } else {
            template[[cov_name]] <- rep(mode_value, curve_grid_length)
          }
        }
      }
      new_design_df <- cbind(new_basis, template)
    } else {
      new_design_df <- new_basis
    }

    new_design <- stats::model.matrix(design_formula, data = new_design_df)
    coef_mat <- stats::coef(fit)
    coef_mat <- coef_mat[curve_genes, , drop = FALSE]
    pred_eta <- coef_mat %*% t(new_design)
    pred_lcpm <- (pred_eta + log(1e6)) / log(2)

    curve_df <- data.frame(
      gene = rep(rownames(pred_lcpm), each = ncol(pred_lcpm)),
      time = rep(time_grid, times = nrow(pred_lcpm)),
      logCPM = as.vector(pred_lcpm),
      condition = label,
      stringsAsFactors = FALSE
    )
    curve_df$status <- ifelse(curve_df$gene %in% tab$gene[tab$significant], "Sig", "NotSig")
  }

  list(
    dge = dge,
    design = design,
    fit = fit,
    test = qlf,
    table = tab,
    curve_data = curve_df,
    time_grid = if (exists("time_grid")) time_grid else numeric(0),
    effective_df = effective_df
  )
}


#' Run pseudobulk spline or time-series differential expression with edgeR
#'
#' Fit spline-based edgeR models to identify genes that change over time within
#' each condition, then optionally return fitted expression curves for
#' visualization.
#'
#' @param pb A pseudobulk `SingleCellExperiment`.
#' @param time_col Sample metadata column containing numeric or ordered time.
#' @param condition_col Optional sample metadata column defining conditions to be
#'   analysed separately, for example treatment.
#' @param covariates Optional additive covariates to include in each
#'   condition-specific model.
#' @param assay_name Assay containing pseudobulk counts.
#' @param normalization Library-size normalization method. Choices are:
#' - `"TMM"`: edgeR TMM normalization.
#' - `"upperquartile"`: edgeR upper-quartile normalization.
#' - `"none"`: no normalization factors beyond raw library sizes.
#' @param df Degrees of freedom for the natural spline basis.
#' @param filter_by_expr Logical indicating whether to run
#'   `edgeR::filterByExpr()`.
#' @param robust Logical passed to `edgeR::glmQLFit()`.
#' @param abundance_trend Logical passed to `edgeR::glmQLFit()`.
#' @param fdr_cutoff FDR cutoff used to define significant genes.
#' @param return_curve_for Which genes should receive fitted spline curves.
#'   Choices are:
#' - `"significant"`: return curves only for significant genes.
#' - `"all"`: return curves for all tested genes.
#' @param curve_grid_length Number of points used when drawing fitted spline
#'   curves.
#'
#' @return A list containing condition-specific edgeR fits, tables of temporal
#'   DE genes, and fitted curve data for plotting.
#' @export
run_edger_spline_de <- function(
    pb,
    time_col = "timepoint",
    condition_col = NULL,
    covariates = NULL,
    assay_name = "counts",
    normalization = c("TMM", "upperquartile", "none"),
    df = 3,
    filter_by_expr = TRUE,
    robust = TRUE,
    abundance_trend = FALSE,
    fdr_cutoff = 0.05,
    return_curve_for = c("significant", "all"),
    curve_grid_length = 200) {
  normalization <- match.arg(normalization)
  return_curve_for <- match.arg(return_curve_for)

  if (!methods::is(pb, "SingleCellExperiment")) {
    stop("`pb` must be a SingleCellExperiment.", call. = FALSE)
  }
  if (!assay_name %in% SummarizedExperiment::assayNames(pb)) {
    stop("Assay `", assay_name, "` is not present in `pb`.", call. = FALSE)
  }

  meta <- as.data.frame(SummarizedExperiment::colData(pb))
  needed_cols <- unique(c(time_col, condition_col, covariates))
  if (!all(needed_cols %in% colnames(meta))) {
    stop("All requested time, condition, and covariate columns must exist in `colData(pb)`.", call. = FALSE)
  }

  keep <- stats::complete.cases(meta[, needed_cols, drop = FALSE])
  if (!any(keep)) {
    stop("No pseudobulk samples have complete values for the requested model.", call. = FALSE)
  }

  meta <- meta[keep, , drop = FALSE]
  counts <- SummarizedExperiment::assay(pb, assay_name)[, keep, drop = FALSE]

  if (is.null(condition_col)) {
    res <- .fit_one_spline_de(
      counts = counts,
      meta = meta,
      time_col = time_col,
      covariates = covariates,
      normalization = normalization,
      df = df,
      filter_by_expr = filter_by_expr,
      robust = robust,
      abundance_trend = abundance_trend,
      fdr_cutoff = fdr_cutoff,
      return_curve_for = return_curve_for,
      curve_grid_length = curve_grid_length,
      label = "all"
    )
    return(list(
      by_condition = list(all = res),
      combined_table = res$table,
      curve_data = res$curve_data,
      settings = list(
        time_col = time_col,
        condition_col = condition_col,
        covariates = covariates,
        assay_name = assay_name,
        normalization = normalization,
        df = df,
        filter_by_expr = filter_by_expr,
        robust = robust,
        abundance_trend = abundance_trend,
        fdr_cutoff = fdr_cutoff,
        return_curve_for = return_curve_for,
        curve_grid_length = curve_grid_length
      )
    ))
  }

  condition_values <- unique(as.character(meta[[condition_col]]))
  condition_values <- condition_values[order(condition_values)]
  res_list <- vector("list", length(condition_values))
  names(res_list) <- condition_values

  for (i in seq_along(condition_values)) {
    cond <- condition_values[i]
    idx <- which(as.character(meta[[condition_col]]) == cond)
    if (length(idx) < max(df + 1L, 4L)) {
      next
    }
    res_list[[i]] <- .fit_one_spline_de(
      counts = counts[, idx, drop = FALSE],
      meta = meta[idx, , drop = FALSE],
      time_col = time_col,
      covariates = covariates,
      normalization = normalization,
      df = df,
      filter_by_expr = filter_by_expr,
      robust = robust,
      abundance_trend = abundance_trend,
      fdr_cutoff = fdr_cutoff,
      return_curve_for = return_curve_for,
      curve_grid_length = curve_grid_length,
      label = cond
    )
  }

  res_list <- Filter(Negate(is.null), res_list)
  if (!length(res_list)) {
    stop("No condition had enough pseudobulk samples for spline DE.", call. = FALSE)
  }

  combined_table <- do.call(rbind, lapply(res_list, `[[`, "table"))
  curve_data <- do.call(rbind, lapply(res_list, `[[`, "curve_data"))

  list(
    by_condition = res_list,
    combined_table = combined_table,
    curve_data = curve_data,
    settings = list(
      time_col = time_col,
      condition_col = condition_col,
      covariates = covariates,
      assay_name = assay_name,
      normalization = normalization,
      df = df,
      filter_by_expr = filter_by_expr,
      robust = robust,
      abundance_trend = abundance_trend,
      fdr_cutoff = fdr_cutoff,
      return_curve_for = return_curve_for,
      curve_grid_length = curve_grid_length
    )
  )
}


#' Plot fitted spline expression curves from time-series DE
#'
#' Plot fitted spline curves for selected genes, optionally overlaid across
#' conditions.
#'
#' @param spline_result Result list returned by `run_edger_spline_de()`.
#' @param genes Optional character vector of genes to plot. If `NULL`, the top
#'   genes ranked by FDR are used.
#' @param top_n Number of genes to plot when `genes = NULL`.
#' @param scales Facet scale setting passed to `ggplot2::facet_wrap()`.
#' @param ncol Optional number of columns used by `ggplot2::facet_wrap()`.
#'
#' @return A `ggplot2` object.
#' @export
plot_time_series_deg_curves <- function(
    spline_result,
    genes = NULL,
    top_n = 12L,
    scales = "free_y",
    ncol = NULL) {
  if (is.null(spline_result$curve_data) || !nrow(spline_result$curve_data)) {
    stop("`spline_result` does not contain any curve data to plot.", call. = FALSE)
  }

  curve_df <- spline_result$curve_data
  table_df <- spline_result$combined_table

  if (is.null(genes)) {
    rank_df <- stats::aggregate(FDR ~ gene, data = table_df, FUN = min)
    rank_df <- rank_df[order(rank_df$FDR, rank_df$gene), , drop = FALSE]
    genes <- utils::head(rank_df$gene, as.integer(top_n))
  }

  curve_df <- curve_df[curve_df$gene %in% genes, , drop = FALSE]
  curve_df$gene <- factor(curve_df$gene, levels = genes)
  curve_df$status <- factor(curve_df$status, levels = c("Sig", "NotSig"))

  colour_label <- if (is.null(spline_result$settings$condition_col)) {
    "condition"
  } else {
    spline_result$settings$condition_col
  }

  ggplot2::ggplot(
    curve_df,
    ggplot2::aes(x = time, y = logCPM, colour = condition, linetype = status, group = interaction(condition, status))
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~ gene, scales = scales, ncol = ncol) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = spline_result$settings$time_col, y = "fitted logCPM", colour = colour_label)
}


#' Run a full SCProkaR differential expression workflow
#'
#' Aggregate single-cell counts into pseudobulk samples and then dispatch to a
#' pairwise or time-series edgeR workflow.
#'
#' @param sce A `SingleCellExperiment`.
#' @param mode Workflow mode. Choices are:
#' - `"pairwise"`: pseudobulk followed by pairwise edgeR contrasts.
#' - `"time_series"`: pseudobulk followed by spline-based time-series edgeR.
#' @param sample_cols Metadata columns used to define pseudobulk samples.
#' @param assay_name Assay containing input counts.
#' @param aggregation Pseudobulk aggregation method. Choices are:
#' - `"sum"`: sum counts within each pseudobulk sample.
#' - `"mean"`: average counts within each pseudobulk sample.
#' @param min_cells Minimum number of cells required per pseudobulk sample.
#' @param normalization Pseudobulk normalization method. Choices are:
#' - `"TMM"`: edgeR TMM normalization.
#' - `"upperquartile"`: edgeR upper-quartile normalization.
#' - `"none"`: no normalization factors beyond raw library sizes.
#' @param ... Additional arguments passed to `run_edger_pairwise_de()` or
#'   `run_edger_spline_de()`, depending on `mode`.
#'
#' @return A list with the pseudobulk object, the DE result, and workflow
#'   settings.
#' @export
run_scproka_de <- function(
    sce,
    mode = c("pairwise", "time_series"),
    sample_cols,
    assay_name = "counts",
    aggregation = c("sum", "mean"),
    min_cells = 10L,
    normalization = c("TMM", "upperquartile", "none"),
    ...) {
  mode <- match.arg(mode)
  aggregation <- match.arg(aggregation)
  normalization <- match.arg(normalization)

  pb <- aggregate_pseudobulk(
    sce = sce,
    sample_cols = sample_cols,
    assay_name = assay_name,
    aggregation = aggregation,
    min_cells = min_cells,
    add_logcounts = TRUE
  )

  pb <- normalize_pseudobulk(
    pb,
    assay_name = "counts",
    method = normalization,
    output_assay = "logcounts"
  )

  de_result <- if (mode == "pairwise") {
    run_edger_pairwise_de(pb, assay_name = "counts", normalization = normalization, ...)
  } else {
    run_edger_spline_de(pb, assay_name = "counts", normalization = normalization, ...)
  }

  list(
    pseudobulk = pb,
    de_result = de_result,
    settings = list(
      mode = mode,
      sample_cols = sample_cols,
      assay_name = assay_name,
      aggregation = aggregation,
      min_cells = min_cells,
      normalization = normalization
    )
  )
}
