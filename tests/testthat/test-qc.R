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
    expect_true(!is.null(S4Vectors::metadata(filtered)$scProkaR$filtering))
})

test_that("RunBacQC handles a single-cell object", {
    ## scuttle::perCellQCMetrics() errors on objects with fewer than two
    ## cells, so RunBacQC takes a direct path there. Guard that it still works.
    counts <- Matrix::Matrix(
        matrix(
            c(4, 6, 2, 8, 0, 5),
            ncol = 1,
            dimnames = list(
                c("rrsA", "rrlB", "rplC", "rpsD", "geneX", "geneY"),
                "cell1"
            )
        ),
        sparse = TRUE
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts)
    )

    qc <- RunBacQC(sce, store = FALSE)

    expect_equal(nrow(qc), 1L)
    expect_equal(qc$total_counts, 25)
    expect_equal(qc$detected_features, 5)
    expect_equal(qc$rrna_counts, 10)
    expect_equal(qc$rrna_fraction, 10 / 25)
    expect_equal(qc$pct_rrna, 40)
    expect_equal(qc$ribo_counts, 10)
})

test_that("RunBacQC reports zero shares for an empty cell", {
    ## cell2 is genuinely empty: every gene is zero.
    counts <- Matrix::Matrix(
        matrix(
            c(4, 6, 2, 0, 0, 0),
            ncol = 2,
            dimnames = list(c("rrsA", "rrlB", "geneX"), c("cell1", "cell2"))
        ),
        sparse = TRUE
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts)
    )

    qc <- RunBacQC(sce, store = FALSE)

    ## An empty library has no meaningful share, so it is reported as a true
    ## zero rather than the NaN that dividing by zero would give.
    expect_equal(qc$total_counts, c(12, 0))
    expect_equal(qc$detected_features, c(3, 0))
    expect_equal(qc$rrna_fraction, c(10 / 12, 0))
    expect_false(any(is.na(qc$rrna_fraction)))
    expect_false(any(is.nan(qc$pct_rrna)))
})
