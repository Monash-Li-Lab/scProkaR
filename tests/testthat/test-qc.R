test_that("RunBacQC computes rRNA fractions", {
    sce <- make_toy_sce()
    sce <- RunBacQC(sce)
    cd <- as.data.frame(SummarizedExperiment::colData(sce))

    expect_true(
        all(c("rrna_fraction", "pct_rrna", "ribo_fraction") %in% colnames(cd))
    )
    expect_equal(cd["cell1", "rrna_counts"], 5)
    expect_equal(cd["cell1", "rrna_fraction"], 5 / cd["cell1", "total_counts"])
})

test_that("RunBacQC can use explicit gene classes", {
    sce <- make_toy_sce()
    qc <- RunBacQC(sce, gene_class_col = "feature_class", store = FALSE)

    expect_s3_class(qc, "data.frame")
    expect_equal(qc["cell2", "rrna_counts"], 0)
    expect_equal(qc["cell2", "ribo_counts"], 8)
})

test_that("FilterBacCells records filtering summary", {
    sce <- make_toy_sce()
    sce <- RunBacQC(sce)
    filtered <- FilterBacCells(sce, min_counts = 3, max_rrna_fraction = 0.95)

    expect_lte(ncol(filtered), ncol(sce))
    expect_true(!is.null(S4Vectors::metadata(filtered)$SCProkaR$filtering))
})
