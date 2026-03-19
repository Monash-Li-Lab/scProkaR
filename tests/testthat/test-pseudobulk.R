test_that("AggregatePseudobulk creates grouped count matrices", {
  sce <- make_toy_sce()
  pb <- AggregatePseudobulk(sce, sample_col = "sample_id", group_cols = c("condition"))

  expect_s3_class(pb, "SCProkaRPseudobulk")
  expect_equal(nrow(pb$counts), nrow(sce))
  expect_true(all(c("sample_id", "condition", "n_cells") %in% colnames(pb$samples)))
})

test_that("RunPseudobulkDE validates contrast inputs", {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    skip("edgeR not installed")
  }

  sce <- make_toy_sce()
  pb <- AggregatePseudobulk(sce, sample_col = "sample_id", group_cols = c("condition"))
  expect_error(
    RunPseudobulkDE(pb, design = ~ condition, contrast = "missing_coef"),
    "contrast"
  )
})

test_that("RunTimeSeriesDE requires time or pseudotime information", {
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    skip("tradeSeq not installed")
  }

  sce <- make_toy_sce()
  expect_error(
    RunTimeSeriesDE(sce, time_col = "missing_time"),
    "time_col"
  )
})
