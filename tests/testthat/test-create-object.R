test_that("CreateBacObject standardizes matrix input", {
    sce <- CreateBacObject(
        toy_counts(),
        cell_metadata = toy_metadata(),
        feature_metadata = toy_feature_metadata(),
        sample_col = "sample_id",
        batch_col = "batch"
    )

    expect_s4_class(sce, "SingleCellExperiment")
    expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
    expect_equal(S4Vectors::metadata(sce)$scProkaR$package, "scProkaR")
    expect_equal(rownames(SummarizedExperiment::colData(sce)), colnames(sce))
})

test_that(".scprokar_as_dgC handles triplet sparse matrices", {
    triplet <- methods::as(
        Matrix::Matrix(toy_counts(), sparse = TRUE),
        "TsparseMatrix"
    )
    out <- .scprokar_as_dgC(triplet)

    expect_s4_class(out, "dgCMatrix")
    expect_identical(dim(out), dim(triplet))
})

test_that("CreateBacObject can run unintegrated preprocessing", {
    skip_if_not_installed("igraph")

    sce <- CreateBacObject(
        toy_counts(),
        cell_metadata = toy_metadata(),
        feature_metadata = toy_feature_metadata(),
        sample_col = "sample_id",
        batch_col = "batch",
        run_unintegrated = TRUE,
        unintegrated_dims = 1:2,
        unintegrated_k = 2
    )

    expect_true("logcounts" %in% SummarizedExperiment::assayNames(sce))
    expect_true("PCA" %in% SingleCellExperiment::reducedDimNames(sce))
    expect_true(
        "unintegrated_clusters" %in%
            colnames(SummarizedExperiment::colData(sce))
    )
    if (requireNamespace("uwot", quietly = TRUE)) {
        expect_true(
            "umap.unintegrated" %in% SingleCellExperiment::reducedDimNames(sce)
        )
    }
})

test_that("CreateBacObject accepts SingleCellExperiment input", {
    skip_if_not_installed("SingleCellExperiment")
    raw_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(raw_counts = toy_counts())
    )
    sce <- CreateBacObject(raw_sce, counts_assay = "raw_counts")
    expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
})

test_that("CreateBacObject recovers Barcode and Symbol fields from SingleCellExperiment input", {
    skip_if_not_installed("SingleCellExperiment")

    counts <- Matrix::Matrix(unname(toy_counts()), sparse = TRUE)
    rownames(counts) <- NULL
    colnames(counts) <- NULL

    raw_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(
            Symbol = rownames(toy_counts()),
            ID = paste0("gene_", seq_len(nrow(toy_counts())))
        ),
        colData = S4Vectors::DataFrame(
            Barcode = colnames(toy_counts()),
            sample_id = toy_metadata()[, "sample_id"],
            batch = toy_metadata()[, "batch"]
        )
    )

    sce <- CreateBacObject(
        raw_sce,
        sample_col = "sample_id",
        batch_col = "batch"
    )

    expect_identical(colnames(sce), colnames(toy_counts()))
    expect_identical(rownames(sce), rownames(toy_counts()))
    expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
})

test_that("CreateBacObject promotes Symbol names so RunBacQC can detect rRNA features", {
    skip_if_not_installed("SingleCellExperiment")

    counts <- Matrix::Matrix(
        matrix(
            c(
                10, 0, 3, 5,
                0, 2, 1, 1,
                4, 1, 0, 0
            ),
            nrow = 3,
            byrow = TRUE
        ),
        sparse = TRUE
    )
    rownames(counts) <- c("gene_id_1", "gene_id_2", "gene_id_3")
    colnames(counts) <- c("cell1", "cell2", "cell3", "cell4")

    raw_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(
            Symbol = c("rrsA", "rplB", "geneX"),
            ID = rownames(counts)
        ),
        colData = S4Vectors::DataFrame(
            Barcode = colnames(counts),
            sample_id = rep("sample1", ncol(counts)),
            batch = rep("batch1", ncol(counts))
        )
    )

    sce <- CreateBacObject(
        raw_sce,
        sample_col = "sample_id",
        batch_col = "batch"
    )
    sce <- RunBacQC(sce)

    expect_identical(rownames(sce), c("rrsA", "rplB", "geneX"))
    expect_true(any(sce$rrna_fraction > 0))
    expect_true(any(sce$ribo_fraction > 0))
})

test_that("CreateBacObject accepts a 10x directory path", {
    tmpdir <- tempfile("tenx_")
    dir.create(tmpdir)
    matrix_dir <- file.path(tmpdir, "filtered_feature_bc_matrix")
    dir.create(matrix_dir)

    counts <- Matrix::Matrix(toy_counts(), sparse = TRUE)
    Matrix::writeMM(counts, file.path(matrix_dir, "matrix.mtx"))
    writeLines(colnames(toy_counts()), file.path(matrix_dir, "barcodes.tsv"))
    feature_table <- cbind(
        rownames(toy_counts()),
        rownames(toy_counts()),
        rep("Gene Expression", nrow(toy_counts()))
    )
    utils::write.table(
        feature_table,
        file = file.path(matrix_dir, "features.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    sce <- CreateBacObject(tmpdir)

    expect_s4_class(sce, "SingleCellExperiment")
    expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
    expect_equal(nrow(sce), nrow(toy_counts()))
    expect_equal(ncol(sce), ncol(toy_counts()))
    expect_equal(
        S4Vectors::metadata(sce)$scProkaR$tenx$source_dir,
        normalizePath(tmpdir, winslash = "/", mustWork = TRUE)
    )
})

test_that("CreateBacObject prefers Symbol-style feature names for direct 10x input", {
    tmpdir <- tempfile("tenx_")
    dir.create(tmpdir)
    matrix_dir <- file.path(tmpdir, "filtered_feature_bc_matrix")
    dir.create(matrix_dir)

    counts <- Matrix::Matrix(toy_counts(), sparse = TRUE)
    Matrix::writeMM(counts, file.path(matrix_dir, "matrix.mtx"))
    writeLines(colnames(toy_counts()), file.path(matrix_dir, "barcodes.tsv"))
    feature_table <- cbind(
        paste0("gene_id_", seq_len(nrow(toy_counts()))),
        rownames(toy_counts()),
        rep("Gene Expression", nrow(toy_counts()))
    )
    utils::write.table(
        feature_table,
        file = file.path(matrix_dir, "features.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    sce <- CreateBacObject(tmpdir)

    expect_identical(rownames(sce), rownames(toy_counts()))
    expect_true("feature_id" %in% colnames(SummarizedExperiment::rowData(sce)))
})

test_that("CreateBacObject can assign fixed sample and batch values for 10x input", {
    tmpdir <- tempfile("tenx_")
    dir.create(tmpdir)
    matrix_dir <- file.path(tmpdir, "filtered_feature_bc_matrix")
    dir.create(matrix_dir)

    counts <- Matrix::Matrix(toy_counts(), sparse = TRUE)
    Matrix::writeMM(counts, file.path(matrix_dir, "matrix.mtx"))
    writeLines(colnames(toy_counts()), file.path(matrix_dir, "barcodes.tsv"))
    feature_table <- cbind(
        rownames(toy_counts()),
        rownames(toy_counts()),
        rep("Gene Expression", nrow(toy_counts()))
    )
    utils::write.table(
        feature_table,
        file = file.path(matrix_dir, "features.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    sce <- CreateBacObject(
        tmpdir,
        sample_id_value = "sample_a",
        batch_value = "batch_a"
    )

    expect_true(all(as.character(sce$sample_id) == "sample_a"))
    expect_true(all(as.character(sce$batch) == "batch_a"))
    expect_identical(
        S4Vectors::metadata(sce)$scProkaR$columns$sample_col,
        "sample_id"
    )
    expect_identical(
        S4Vectors::metadata(sce)$scProkaR$columns$batch_col,
        "batch"
    )
})

test_that("10x import helper can recover cell names from Barcode metadata", {
    ## Build an object whose cells are unnamed but whose colData carries the
    ## original barcodes, which is what a raw 10x import looks like.
    unnamed_counts <- Matrix::Matrix(toy_counts(), sparse = TRUE)
    colnames(unnamed_counts) <- NULL

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = unnamed_counts),
        colData = S4Vectors::DataFrame(Barcode = colnames(toy_counts()))
    )
    colnames(sce) <- NULL

    sce <- .scprokar_ensure_cell_names(sce)

    expect_identical(colnames(sce), colnames(toy_counts()))
})

test_that("CreateBacObject accepts Seurat input when SeuratObject is available", {
    skip_if_not_installed("SeuratObject")

    seu <- SeuratObject::CreateSeuratObject(counts = toy_counts(), assay = "RNA")
    seu$sample_id <- toy_metadata()[colnames(seu), "sample_id"]
    seu$batch <- toy_metadata()[colnames(seu), "batch"]

    pca_emb <- matrix(
        seq_len(ncol(seu) * 2),
        ncol = 2,
        dimnames = list(colnames(seu), c("PC_1", "PC_2"))
    )
    umap_emb <- matrix(
        seq_len(ncol(seu) * 2),
        ncol = 2,
        dimnames = list(colnames(seu), c("UMAP_1", "UMAP_2"))
    )
    seu[["pca"]] <- SeuratObject::CreateDimReducObject(
        embeddings = pca_emb,
        assay = "RNA",
        key = "PC_"
    )
    seu[["umap"]] <- SeuratObject::CreateDimReducObject(
        embeddings = umap_emb,
        assay = "RNA",
        key = "UMAP_"
    )

    sce <- CreateBacObject(seu, seurat_assay = "RNA", sample_col = "sample_id", batch_col = "batch")

    expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
    expect_true(
        all(c("pca", "umap") %in% SingleCellExperiment::reducedDimNames(sce))
    )
    expect_equal(S4Vectors::metadata(sce)$scProkaR$seurat$assay, "RNA")
})

test_that("MergeBacObjects merges per-sample objects safely", {
    counts1 <- toy_counts()[, 1:3, drop = FALSE]
    counts2 <- toy_counts()[, 4:6, drop = FALSE]
    colnames(counts2) <- colnames(counts1)

    meta1 <- toy_metadata()[1:3, , drop = FALSE]
    meta2 <- toy_metadata()[4:6, , drop = FALSE]
    rownames(meta2) <- colnames(counts2)

    sce1 <- CreateBacObject(
        counts1,
        cell_metadata = meta1,
        feature_metadata = toy_feature_metadata(),
        sample_col = "sample_id",
        batch_col = "batch"
    )
    SingleCellExperiment::reducedDim(sce1, "pca") <- matrix(
        rnorm(ncol(sce1) * 2),
        ncol = 2,
        dimnames = list(colnames(sce1), c("PC_1", "PC_2"))
    )
    SingleCellExperiment::reducedDim(sce1, "umap") <- matrix(
        rnorm(ncol(sce1) * 2),
        ncol = 2,
        dimnames = list(colnames(sce1), c("UMAP_1", "UMAP_2"))
    )
    sce2 <- CreateBacObject(
        counts2,
        cell_metadata = meta2,
        feature_metadata = toy_feature_metadata(),
        sample_col = "sample_id",
        batch_col = "batch"
    )
    SingleCellExperiment::reducedDim(sce2, "pca") <- matrix(
        rnorm(ncol(sce2) * 2),
        ncol = 2,
        dimnames = list(colnames(sce2), c("PC_1", "PC_2"))
    )
    SingleCellExperiment::reducedDim(sce2, "umap") <- matrix(
        rnorm(ncol(sce2) * 2),
        ncol = 2,
        dimnames = list(colnames(sce2), c("UMAP_1", "UMAP_2"))
    )

    merged <- MergeBacObjects(sce1, sce2)

    expect_s4_class(merged, "SingleCellExperiment")
    expect_equal(ncol(merged), ncol(sce1) + ncol(sce2))
    expect_equal(rownames(merged), rownames(sce1))
    expect_identical(anyDuplicated(colnames(merged)), 0L)
    expect_true(
        "original_cell_id" %in% colnames(SummarizedExperiment::colData(merged))
    )
    expect_equal(S4Vectors::metadata(merged)$scProkaR$merge$n_objects, 2)
    expect_true(
        all(c("pca", "umap") %in%
            SingleCellExperiment::reducedDimNames(merged))
    )
    expect_equal(
        nrow(SingleCellExperiment::reducedDim(merged, "pca")),
        ncol(merged)
    )
    expect_equal(
        nrow(SingleCellExperiment::reducedDim(merged, "umap")),
        ncol(merged)
    )
    expect_equal(
        ncol(SingleCellExperiment::reducedDim(merged, "pca")),
        ncol(SingleCellExperiment::reducedDim(sce1, "pca"))
    )
})

test_that("MergeBacObjects merges lists of more than two objects", {
    sce1 <- CreateBacObject(
        toy_counts()[, 1:2, drop = FALSE],
        cell_metadata = toy_metadata()[1:2, , drop = FALSE],
        feature_metadata = toy_feature_metadata(),
        sample_col = "sample_id",
        batch_col = "batch"
    )
    sce2 <- CreateBacObject(
        toy_counts()[, 3:4, drop = FALSE],
        cell_metadata = toy_metadata()[3:4, , drop = FALSE],
        feature_metadata = toy_feature_metadata(),
        sample_col = "sample_id",
        batch_col = "batch"
    )
    sce3 <- CreateBacObject(
        toy_counts()[, 5:6, drop = FALSE],
        cell_metadata = toy_metadata()[5:6, , drop = FALSE],
        feature_metadata = toy_feature_metadata(),
        sample_col = "sample_id",
        batch_col = "batch"
    )

    merged <- MergeBacObjects(objects = list(sce1, sce2, sce3), gene_mode = "intersect")

    expect_s4_class(merged, "SingleCellExperiment")
    expect_equal(ncol(merged), 6)
    expect_equal(nrow(SummarizedExperiment::colData(merged)), 6)
    expect_equal(ncol(SummarizedExperiment::assay(merged, "counts")), 6)
})
