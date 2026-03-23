test_that("BenchmarkIntegration returns scores and rankings", {
  sce <- make_toy_sce()
  sce <- RunBacQC(sce)
  sce <- SCProkaR:::.scprokar_prepare_reduction(sce, dims = 1:2)
  SingleCellExperiment::reducedDim(sce, "integrated_mock") <- SingleCellExperiment::reducedDim(sce, "PCA")

  result <- BenchmarkIntegration(
    sce,
    batch_col = "batch",
    label_col = "cell_type",
    methods = "integrated_mock"
  )

  expect_true(all(c("scores", "ranking") %in% names(result)))
  expect_true(nrow(result$scores) > 0)
  expect_true(nrow(result$ranking) >= 1)
  expect_true(all(c(
    "overall",
    "key_metrics",
    "batch_removal",
    "bio_conservation",
    "tradeoff",
    "heatmap"
  ) %in% names(result$plots)))
  expect_s3_class(result$plots$key_metrics, "ggplot")
  expect_s3_class(result$plots$tradeoff, "ggplot")
})
