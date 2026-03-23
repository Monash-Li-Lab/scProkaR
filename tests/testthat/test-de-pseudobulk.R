test_that("aggregate_pseudobulk and normalize_pseudobulk create sample-level objects", {
  sce <- make_toy_multidrug_sce(n_cells = 720, n_features = 60, seed = 11)
  sce <- assign_toy_pseudoreplicates(sce, n_reps = 3, seed = 11)

  pb <- aggregate_pseudobulk(
    sce,
    sample_cols = c("pb_sample_id", "condition", "timepoint", "pb_replicate"),
    aggregation = "sum",
    min_cells = 10
  )

  expect_s4_class(pb, "SingleCellExperiment")
  expect_true("counts" %in% SummarizedExperiment::assayNames(pb))
  expect_true(all(c("condition", "timepoint", "pb_replicate", "ncells", "lib.size") %in%
    colnames(SummarizedExperiment::colData(pb))))
  expect_true(all(SummarizedExperiment::colData(pb)$ncells >= 10))

  pb <- normalize_pseudobulk(pb, method = "TMM")
  expect_true("logcounts" %in% SummarizedExperiment::assayNames(pb))
  expect_true(all(c("lib.size", "norm.factors") %in% colnames(SummarizedExperiment::colData(pb))))
})

test_that("filter_pseudobulk_samples removes small pseudobulk samples", {
  pb <- make_toy_pseudobulk(n_cells = 720, n_features = 60, n_reps = 3, seed = 12)
  filtered <- filter_pseudobulk_samples(pb, min_cells = 15, min_lib_size = 1)

  expect_s4_class(filtered, "SingleCellExperiment")
  expect_true(ncol(filtered) <= ncol(pb))
  expect_true(all(SummarizedExperiment::colData(filtered)$ncells >= 15))
})
