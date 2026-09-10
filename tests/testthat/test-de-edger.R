test_that("run_edger_pairwise_de returns contrast tables and plotting helpers work", {
    pb <- make_toy_pseudobulk(n_cells = 720, n_features = 60, n_reps = 3, seed = 21)

    de <- run_edger_pairwise_de(
        pb,
        group_col = "timepoint",
        covariates = "condition",
        contrast_type = "reference",
        reference_level = "0",
        normalization = "TMM"
    )

    expect_true(
        all(c("dge", "design", "fit", "tables", "combined_table") %in%
            names(de))
    )
    expect_true("24_vs_0" %in% names(de$tables))
    expect_true(nrow(de$combined_table) > 0)
    expect_true(
        all(c("gene", "logFC", "FDR", "contrast") %in%
            colnames(de$combined_table))
    )

    volcano <- plot_pairwise_de_volcano(de, contrast = "24_vs_0", top_n_labels = 5)
    ma <- plot_pairwise_de_ma(de, contrast = "24_vs_0")

    expect_s3_class(volcano, "ggplot")
    expect_s3_class(ma, "ggplot")
})

test_that("run_edger_spline_de and run_scproka_de return time-series results", {
    pb <- make_toy_pseudobulk(n_cells = 720, n_features = 60, n_reps = 3, seed = 22)

    spline <- run_edger_spline_de(
        pb,
        time_col = "timepoint",
        condition_col = "condition",
        df = 3,
        normalization = "TMM",
        return_curve_for = "all"
    )

    expect_true(
        all(c("by_condition", "combined_table", "curve_data") %in%
            names(spline))
    )
    expect_true(
        all(c("drug_A", "drug_B", "drug_C") %in% names(spline$by_condition))
    )
    expect_true(nrow(spline$combined_table) > 0)
    expect_true(nrow(spline$curve_data) > 0)

    curve_plot <- plot_time_series_deg_curves(
        spline,
        genes = unique(utils::head(spline$combined_table$gene, 4)),
        ncol = 2
    )
    expect_s3_class(curve_plot, "ggplot")

    sce <- make_toy_multidrug_sce(n_cells = 720, n_features = 60, seed = 23)
    sce <- assign_toy_pseudoreplicates(sce, n_reps = 3, seed = 23)
    wrapped <- run_scproka_de(
        sce,
        mode = "pairwise",
        sample_cols = c("pb_sample_id", "condition", "timepoint", "pb_replicate"),
        min_cells = 10,
        group_col = "timepoint",
        covariates = "condition",
        contrast_type = "reference",
        reference_level = "0"
    )

    expect_true(
        all(c("pseudobulk", "de_result", "settings") %in% names(wrapped))
    )
    expect_s4_class(wrapped$pseudobulk, "SingleCellExperiment")
    expect_true("24_vs_0" %in% names(wrapped$de_result$tables))
})
