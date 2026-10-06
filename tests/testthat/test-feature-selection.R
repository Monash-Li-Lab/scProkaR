test_that("feature selection respects cached HVGs before modelling", {
    counts <- matrix(
        seq_len(20), nrow = 5,
        dimnames = list(paste0("gene", seq_len(5)), paste0("cell", seq_len(4)))
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts)
    )
    SummarizedExperiment::rowData(sce)$is_hvg <-
        c(FALSE, TRUE, NA, FALSE, TRUE)

    expect_identical(
        .scprokar_select_features(sce, nfeatures = 1),
        c("gene2", "gene5")
    )
    expect_identical(
        .scprokar_select_features(sce, feature_set = "all"),
        rownames(sce)
    )
})

test_that("feature selection retains sparse variance fallback", {
    logcounts <- Matrix::Matrix(
        rbind(
            constant = c(2, 2, 2, 2),
            variable = c(0, 0, 0, 8),
            medium = c(0, 2, 0, 2),
            silent = c(0, 0, 0, 0)
        ),
        sparse = TRUE
    )
    colnames(logcounts) <- paste0("cell", seq_len(ncol(logcounts)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = logcounts)
    )

    ## Exercise an installation without the optional backend even when the
    ## check environment has scrapper installed.
    select_without_backend <- .scprokar_select_features
    environment(select_without_backend) <- list2env(
        list(requireNamespace = function(package, quietly) FALSE),
        parent = environment(.scprokar_select_features)
    )

    expect_identical(
        select_without_backend(sce, nfeatures = 2),
        c("variable", "medium")
    )
    expect_length(select_without_backend(sce, nfeatures = 100), nrow(sce))
})

test_that("feature selection ranks scrapper variance residuals", {
    skip_if_not_installed("scrapper", minimum_version = "1.0.0")

    set.seed(1206)
    ngenes <- 300L
    ncells <- 40L
    means <- exp(seq(log(0.1), log(20), length.out = ngenes))
    counts <- matrix(
        stats::rpois(ngenes * ncells, lambda = rep(means, ncells)),
        nrow = ngenes,
        dimnames = list(paste0("gene", seq_len(ngenes)),
                        paste0("cell", seq_len(ncells)))
    )
    logcounts <- Matrix::Matrix(log1p(counts), sparse = TRUE)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = logcounts)
    )

    fitted <- scrapper::modelGeneVariances(logcounts)$statistics
    expected <- rownames(sce)[
        order(fitted$residuals, decreasing = TRUE)[seq_len(20)]
    ]
    selected <- NULL
    expect_warning(
        selected <- .scprokar_select_features(sce, nfeatures = 20), NA
    )
    expect_identical(selected, expected)
    expect_false(identical(
        selected,
        rownames(sce)[order(fitted$variances, decreasing = TRUE)[seq_len(20)]]
    ))
})
