existing_test_embedding <- function(sce, offset = 0) {
    matrix(
        sin(seq_len(ncol(sce) * 3) + offset),
        ncol = 3,
        dimnames = list(colnames(sce), paste0("latent_", seq_len(3)))
    )
}

test_that("existing embeddings are detected and registered with provenance", {
    sce <- make_toy_sce()
    cca <- existing_test_embedding(sce)
    scvi <- existing_test_embedding(sce, offset = 1)
    SingleCellExperiment::reducedDim(sce, "seurat_cca") <- cca
    SingleCellExperiment::reducedDim(sce, "X_scVI") <- scvi
    SingleCellExperiment::reducedDim(sce, "PCA") <- cca
    SingleCellExperiment::reducedDim(sce, "X_umap") <- cca[, seq_len(2)]
    SingleCellExperiment::reducedDim(sce, "umap_integrated_scvi") <-
        scvi[, seq_len(2)]

    out <- RegisterExistingIntegrationEmbeddings(
        sce, metadata = list(scvi = list(version = "test"))
    )
    records <- S4Vectors::metadata(out)$scProkaR$integration$results
    expect_setequal(names(records), c("cca", "scvi"))
    expect_equal(records$cca$original_reduction, "seurat_cca")
    expect_equal(records$scvi$source, "external_existing_reduction")
    expect_equal(records$scvi$version, "test")
    expect_equal(SingleCellExperiment::reducedDim(out, "integrated_cca"), cca)
    expect_equal(SingleCellExperiment::reducedDim(out, "integrated_scvi"), scvi)
    expect_equal(SingleCellExperiment::reducedDim(out, "X_scVI"), scvi)
    expect_null(S4Vectors::metadata(out)$SCProkaR)
})

test_that("existing embeddings support references and custom prefixes", {
    sce <- make_toy_sce()
    embedding <- existing_test_embedding(sce)
    SingleCellExperiment::reducedDim(sce, "X_scVI") <- embedding
    out <- RegisterExistingIntegrationEmbeddings(
        sce, reductions = "X_scVI", method_names = "latent", copy = FALSE,
        metadata = list(software = "external")
    )
    expect_identical(SingleCellExperiment::reducedDimNames(out), "X_scVI")
    record <- S4Vectors::metadata(out)$scProkaR$integration$results$latent
    expect_equal(record$reduction, "X_scVI")
    expect_equal(record$software, "external")

    copied <- RegisterExistingIntegrationEmbeddings(
        sce, reductions = "X_scVI", prefix = "corrected_"
    )
    expect_equal(
        SingleCellExperiment::reducedDim(copied, "corrected_scvi"), embedding
    )
    repeated <- RegisterExistingIntegrationEmbeddings(
        copied, reductions = "corrected_scvi", prefix = "corrected_",
        method_names = "scvi"
    )
    expect_identical(
        SingleCellExperiment::reducedDimNames(repeated),
        SingleCellExperiment::reducedDimNames(copied)
    )
})

test_that("existing embedding registration rejects invalid selections", {
    expect_error(RegisterExistingIntegrationEmbeddings(list()),
        "must be a SingleCellExperiment")
    sce <- make_toy_sce()
    expect_error(RegisterExistingIntegrationEmbeddings(sce),
        "No reduced dimensions")
    SingleCellExperiment::reducedDim(sce, "PCA") <-
        existing_test_embedding(sce)
    expect_error(RegisterExistingIntegrationEmbeddings(sce),
        "No likely external integration reductions")
    expect_error(RegisterExistingIntegrationEmbeddings(sce, "absent"),
        "were not found")
    expect_error(RegisterExistingIntegrationEmbeddings(sce, character()),
        "non-empty names")
    expect_error(RegisterExistingIntegrationEmbeddings(
        sce, "PCA", method_names = c("one", "two")), "same length")
    expect_error(RegisterExistingIntegrationEmbeddings(
        sce, c("PCA", "PCA")), "unique, non-empty")
    expect_error(RegisterExistingIntegrationEmbeddings(
        sce, "PCA", method_names = ""), "unique, non-empty")
})

test_that("registration protects embeddings unless overwrite is explicit", {
    sce <- make_toy_sce()
    embedding <- existing_test_embedding(sce)
    sce <- RegisterIntegrationEmbedding(sce, embedding, "custom")
    changed <- embedding
    changed[1, 1] <- changed[1, 1] + 1e-10
    expect_error(RegisterIntegrationEmbedding(sce, changed, "custom"),
        "overwrite = TRUE", fixed = TRUE)
    expect_equal(
        SingleCellExperiment::reducedDim(sce, "integrated_custom"), embedding
    )

    reordered <- embedding[rev(seq_len(nrow(embedding))), , drop = FALSE]
    same <- RegisterIntegrationEmbedding(sce, reordered, "custom")
    expect_equal(
        SingleCellExperiment::reducedDim(same, "integrated_custom"), embedding
    )
    out <- RegisterIntegrationEmbedding(
        sce, changed, "custom", overwrite = TRUE,
        metadata = list(source = "imported", version = "test")
    )
    expect_identical(
        SingleCellExperiment::reducedDim(out, "integrated_custom"), changed
    )
    record <- S4Vectors::metadata(out)$scProkaR$integration$results$custom
    expect_equal(record$source, "imported")
    expect_equal(sum(names(record) == "source"), 1)
    expect_equal(record$version, "test")
})

test_that("copying existing reductions respects overwrite protection", {
    sce <- make_toy_sce()
    original <- existing_test_embedding(sce)
    replacement <- existing_test_embedding(sce, offset = 1)
    SingleCellExperiment::reducedDim(sce, "integrated_scvi") <- original
    SingleCellExperiment::reducedDim(sce, "X_scVI") <- replacement
    expect_error(RegisterExistingIntegrationEmbeddings(sce, "X_scVI"),
        "already exists")
    out <- RegisterExistingIntegrationEmbeddings(
        sce, "X_scVI", overwrite = TRUE
    )
    expect_equal(
        SingleCellExperiment::reducedDim(out, "integrated_scvi"), replacement
    )
})

test_that("benchmarking includes registered and discovered embeddings", {
    sce <- make_toy_sce()
    embedding <- existing_test_embedding(sce)
    SingleCellExperiment::reducedDim(sce, "seurat_cca") <- embedding
    sce <- RegisterExistingIntegrationEmbeddings(sce, "seurat_cca")
    SingleCellExperiment::reducedDim(sce, "integrated_scanorama") <-
        existing_test_embedding(sce, offset = 1)
    SingleCellExperiment::reducedDim(sce, "PCA") <- embedding

    result <- BenchmarkIntegration(
        sce, batch_col = "batch", label_col = "cell_type", return_plots = FALSE
    )
    expect_setequal(unique(result$scores$method), c("cca", "scanorama"))
    expect_setequal(unique(result$scores$reduction),
        c("integrated_cca", "integrated_scanorama"))
})

test_that("benchmark method and reduction matching ignores case", {
    sce <- make_toy_sce()
    embedding <- existing_test_embedding(sce)
    sce <- RegisterIntegrationEmbedding(sce, embedding, "Harmony")
    SingleCellExperiment::reducedDim(sce, "integrated_scanorama") <- embedding
    SingleCellExperiment::reducedDim(sce, "X_scVI") <- embedding

    result <- BenchmarkIntegration(
        sce, batch_col = "batch", methods = c("hARMONY", "SCANORAMA", "x_ScVi"),
        return_plots = FALSE
    )
    expect_setequal(unique(result$scores$method),
        c("Harmony", "scanorama", "X_scVI"))
    expect_setequal(unique(result$scores$reduction),
        c("integrated_Harmony", "integrated_scanorama", "X_scVI"))

    by_reduction <- BenchmarkIntegration(
        sce, batch_col = "batch", methods = "INTEGRATED_SCANORAMA",
        return_plots = FALSE
    )
    expect_identical(unique(by_reduction$scores$reduction),
        "integrated_scanorama")
})

test_that("discovery preserves registered aliases and fallback reductions", {
    sce <- make_toy_sce()
    embedding <- existing_test_embedding(sce)
    SingleCellExperiment::reducedDim(sce, "PCA") <- embedding
    expect_identical(.scprokar_reduction_lookup(sce), list(PCA = "PCA"))

    sce <- RegisterIntegrationEmbedding(
        sce, embedding, "scanorama", reduction_name = "latent"
    )
    SingleCellExperiment::reducedDim(sce, "integrated_scanorama") <- embedding
    expect_identical(.scprokar_reduction_lookup(sce), list(
        scanorama = "latent", integrated_scanorama = "integrated_scanorama"
    ))
})

test_that("exact reduction names take precedence over case-folded methods", {
    sce <- make_toy_sce()
    embedding <- existing_test_embedding(sce)
    SingleCellExperiment::reducedDim(sce, "PCA") <- embedding
    sce <- RegisterIntegrationEmbedding(sce, embedding, "pca")
    expect_identical(.scprokar_reduction_lookup(sce, "PCA"), list(PCA = "PCA"))
    expect_identical(.scprokar_reduction_lookup(sce, "pca"),
        list(pca = "integrated_pca"))
})

test_that("an exact reducedDim name outranks an auto-discovered alias", {
    set.seed(1)
    sce <- simulate_tata_multidrug_sce(
        n_cells = 300, n_features = 60, n_pcs = 10
    )
    pca <- SingleCellExperiment::reducedDim(sce, "PCA")

    ## A reduction literally named "mnn" alongside a different, unregistered
    ## "integrated_mnn". Asking for "mnn" must score the one named "mnn".
    SingleCellExperiment::reducedDim(sce, "mnn") <- pca[, seq_len(4)]
    SingleCellExperiment::reducedDim(sce, "integrated_mnn") <- pca[, 5:8]

    lookup <- .scprokar_reduction_lookup(sce, methods = "mnn")

    expect_identical(unname(unlist(lookup)), "mnn")
})

test_that("a registered method name still resolves to its reduction", {
    set.seed(1)
    sce <- simulate_tata_multidrug_sce(
        n_cells = 300, n_features = 60, n_pcs = 10
    )
    pca <- SingleCellExperiment::reducedDim(sce, "PCA")
    sce <- RegisterIntegrationEmbedding(
        sce, pca[, seq_len(4)],
        method_name = "mnn"
    )

    lookup <- .scprokar_reduction_lookup(sce, methods = "mnn")

    expect_identical(unname(unlist(lookup)), "integrated_mnn")
})

test_that("auto-detection leaves scProkaR's own provenance intact", {
    set.seed(1)
    sce <- simulate_tata_multidrug_sce(
        n_cells = 300, n_features = 60, n_pcs = 10
    )
    pca <- SingleCellExperiment::reducedDim(sce, "PCA")
    sce <- RegisterIntegrationEmbedding(
        sce, pca[, seq_len(4)],
        method_name = "mnn",
        metadata = list(batch_col = "condition", feature_set = "hvg")
    )
    SingleCellExperiment::reducedDim(sce, "X_scVI") <- pca[, seq_len(3)]

    out <- RegisterExistingIntegrationEmbeddings(sce)
    record <- S4Vectors::metadata(out)$scProkaR$integration$results$mnn

    ## The registry entry written by the earlier registration survives.
    expect_identical(record$batch_col, "condition")
    expect_identical(record$feature_set, "hvg")
    expect_true("scvi" %in% names(
        S4Vectors::metadata(out)$scProkaR$integration$results
    ))
})

test_that("auto-detection is idempotent", {
    set.seed(1)
    sce <- simulate_tata_multidrug_sce(
        n_cells = 300, n_features = 60, n_pcs = 10
    )
    pca <- SingleCellExperiment::reducedDim(sce, "PCA")
    SingleCellExperiment::reducedDim(sce, "X_scVI") <- pca[, seq_len(3)]

    once <- RegisterExistingIntegrationEmbeddings(sce)

    expect_error(RegisterExistingIntegrationEmbeddings(once), NA)
})

test_that("visualisation layouts are not mistaken for latent spaces", {
    detected <- .scprokar_guess_external_integration_reductions(c(
        "X_pca", "X_umap", "X_tsne", "X_diffmap", "X_draw_graph_fa",
        "X_phate", "X_densmap", "harmony_umap",
        "X_scVI", "seurat_cca", "scanorama_latent"
    ))

    expect_setequal(detected, c("X_scVI", "seurat_cca", "scanorama_latent"))
})

test_that("metadata keys beginning with source are not promoted", {
    set.seed(1)
    sce <- simulate_tata_multidrug_sce(
        n_cells = 300, n_features = 60, n_pcs = 10
    )
    pca <- SingleCellExperiment::reducedDim(sce, "PCA")

    out <- RegisterIntegrationEmbedding(
        sce, pca[, seq_len(3)],
        method_name = "m1",
        metadata = list(source_file = "scvi.h5ad")
    )
    record <- S4Vectors::metadata(out)$scProkaR$integration$results$m1

    ## `$` partial matching would have promoted source_file into source.
    expect_identical(record$source, "external")
    expect_identical(record$source_file, "scvi.h5ad")
})
