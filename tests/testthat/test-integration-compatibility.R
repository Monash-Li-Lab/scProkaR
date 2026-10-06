reference_fast_mnn <- function(...) {
    withCallingHandlers(
        batchelor::fastMNN(...),
        warning = function(w) {
            ## The reference intentionally exercises the old dependency path.
            if (inherits(w, "deprecatedWarning") &&
                    grepl("normalizeCounts", conditionMessage(w), fixed = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }
    )
}

test_that("MNN preprocessing preserves the original batch correction", {
    skip_if_not_installed("batchelor")
    skip_if_not_installed("BiocSingular")
    set.seed(103)
    logx <- matrix(log1p(rpois(80 * 54, lambda = 4)), nrow = 80)
    logx[, 1] <- 0
    dimnames(logx) <- list(paste0("g", 1:80), paste0("c", 1:54))
    features <- rownames(logx)[seq(1, 80, by = 3)]

    for (sparse in c(FALSE, TRUE)) {
        x <- if (sparse) Matrix::Matrix(logx, sparse = TRUE) else logx
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(logcounts = x)
        )
        batches <- lapply(split(seq_len(ncol(sce)), rep(1:3, each = 18)),
            function(i) sce[, i])
        for (correct_all in c(FALSE, TRUE)) {
            args <- list(
                k = 4, BSPARAM = BiocSingular::ExactParam(),
                correct.all = correct_all, get.variance = TRUE,
                merge.order = c(3, 1, 2)
            )
            reference <- do.call(reference_fast_mnn,
                c(batches, list(subset.row = features, d = 5), args))
            prepared <- .scprokar_prepare_mnn_batches(
                batches, features = features, args = args
            )
            expect_warning(actual <- do.call(batchelor::fastMNN,
                c(prepared$batches,
                    list(subset.row = features, d = 5), prepared$args)), NA)
            expect_equal(
                SingleCellExperiment::reducedDim(actual, "corrected"),
                SingleCellExperiment::reducedDim(reference, "corrected"),
                tolerance = 1e-10
            )
            expect_equal(
                SummarizedExperiment::assay(actual, "reconstructed"),
                SummarizedExperiment::assay(reference, "reconstructed"),
                tolerance = 1e-10
            )
            expect_equal(S4Vectors::metadata(actual),
                S4Vectors::metadata(reference), tolerance = 1e-10)
        }
    }
})

test_that("MNN preprocessing respects an alternative assay and opt-out", {
    skip_if_not_installed("batchelor")
    x <- Matrix::Matrix(matrix(c(0, 0, 3, 4, 6, 8), nrow = 2),
        sparse = TRUE)
    dimnames(x) <- list(c("g1", "g2"), c("c1", "c2", "c3"))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = x * 2, alternative = x)
    )
    args <- list(assay.type = "alternative")
    prepared <- .scprokar_prepare_mnn_batches(list(sce), "g2", args)
    normalized <- SummarizedExperiment::assay(
        prepared$batches[[1]], "alternative"
    )
    expect_s4_class(normalized, "dgCMatrix")
    expect_equal(as.matrix(normalized),
        matrix(c(0, 0, 0.75, 1, 0.75, 1), nrow = 2, dimnames = dimnames(x)))
    expect_equal(SummarizedExperiment::assay(
        prepared$batches[[1]], "logcounts"), x * 2)
    expect_false(prepared$args$cos.norm)

    args$cos.norm <- FALSE
    unchanged <- .scprokar_prepare_mnn_batches(list(sce), "g2", args)
    expect_identical(unchanged$batches, list(sce))
    expect_identical(unchanged$args, args)
})

test_that("Harmony legacy fallback resolves only an available export", {
    runner <- .scprokar_run_harmony
    shim <- new.env(parent = environment(runner))
    environment(runner) <- shim
    shim$getNamespaceExports <- function(ns) "HarmonyMatrix"
    resolved <- FALSE
    shim$getExportedValue <- function(ns, name) {
        expect_identical(ns, "harmony")
        expect_identical(name, "HarmonyMatrix")
        resolved <<- TRUE
        function(data_mat, meta_data, vars_use, do_pca, nclust, verbose) {
            expect_false(do_pca)
            data_mat + 1
        }
    }
    x <- matrix(seq_len(12), nrow = 6)
    expect_equal(runner(x, data.frame(batch = rep(1:2, each = 3)),
        "batch"), x + 1)
    expect_true(resolved)

    shim$getNamespaceExports <- function(ns) character()
    shim$getExportedValue <- function(...) {
        stop("An unavailable export should never be resolved")
    }
    expect_error(runner(x, data.frame(batch = rep(1:2, each = 3)),
        "batch"), "Harmony integration failed")
})
