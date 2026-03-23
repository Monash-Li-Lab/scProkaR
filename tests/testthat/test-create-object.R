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
  expect_equal(S4Vectors::metadata(sce)$SCProkaR$package, "SCProkaR")
  expect_equal(rownames(SummarizedExperiment::colData(sce)), colnames(sce))
})

test_that("CreateBacObject accepts SingleCellExperiment input", {
  skip_if_not_installed("SingleCellExperiment")
  raw_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(raw_counts = toy_counts())
  )
  sce <- CreateBacObject(raw_sce, counts_assay = "raw_counts")
  expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
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
  expect_true(all(c("pca", "umap") %in% SingleCellExperiment::reducedDimNames(sce)))
  expect_equal(S4Vectors::metadata(sce)$SCProkaR$seurat$assay, "RNA")
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
  expect_true("original_cell_id" %in% colnames(SummarizedExperiment::colData(merged)))
  expect_equal(S4Vectors::metadata(merged)$SCProkaR$merge$n_objects, 2)
  expect_true(all(c("pca", "umap") %in% SingleCellExperiment::reducedDimNames(merged)))
  expect_equal(nrow(SingleCellExperiment::reducedDim(merged, "pca")), ncol(merged))
  expect_equal(nrow(SingleCellExperiment::reducedDim(merged, "umap")), ncol(merged))
  expect_equal(ncol(SingleCellExperiment::reducedDim(merged, "pca")), ncol(SingleCellExperiment::reducedDim(sce1, "pca")))
})
