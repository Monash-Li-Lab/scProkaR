make_trajectory_sce <- function() {
  set.seed(42)
  doses <- c(0, 0.5, 2)
  times <- c(0, 1, 4, 7)
  reps <- 8
  design <- expand.grid(
    dose = doses,
    time = times,
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  n_cells <- nrow(design)
  cells <- paste0("traj_cell_", seq_len(n_cells))
  cluster_id <- ifelse(design$time <= 1, "c1", ifelse(design$dose == 0, "c2", "c3"))
  cluster_centers <- rbind(c(-2, -1), c(0.2, 1.4), c(2.1, 0.4))
  rownames(cluster_centers) <- c("c1", "c2", "c3")

  emb <- matrix(0, nrow = n_cells, ncol = 2, dimnames = list(cells, c("FDL_1", "FDL_2")))
  for (i in seq_len(n_cells)) {
    center <- cluster_centers[cluster_id[i], ]
    emb[i, 1] <- center[1] + design$time[i] * 0.25 + design$dose[i] * 0.5 + rnorm(1, sd = 0.12)
    emb[i, 2] <- center[2] + design$time[i] * 0.12 - design$dose[i] * 0.25 + rnorm(1, sd = 0.12)
  }

  genes <- paste0("gene", seq_len(18))
  counts <- matrix(rpois(length(genes) * n_cells, lambda = 4), nrow = length(genes))
  rownames(counts) <- genes
  colnames(counts) <- cells

  meta <- data.frame(
    seurat_clusters = cluster_id,
    time = design$time,
    dose = design$dose,
    trajectory_state = "real",
    row.names = cells,
    stringsAsFactors = FALSE
  )
  meta$trajectory_state[seq_len(6)] <- "synthetic"

  sce <- CreateBacObject(counts, cell_metadata = meta)
  SingleCellExperiment::reducedDim(sce, "fdl") <- emb
  sce
}

test_that("main trajectory runner stores paga, density, and composition", {
  sce <- make_trajectory_sce()
  out <- RunTrajectory(
    sce,
    reduction = c("fdl", "X_fdl"),
    cluster_col = "seurat_clusters",
    time_col = "time",
    dose_col = "dose",
    trajectory_col = "trajectory_state"
  )

  traj <- S4Vectors::metadata(out)$SCProkaR$trajectory
  expect_true(!is.null(traj$paga$nodes))
  expect_true(!is.null(traj$paga$edges))
  expect_true(!is.null(traj$composition))
  expect_true(!is.null(traj$density$points))
})

test_that("PAGA nodes have one row per cluster and weights are bounded", {
  sce <- make_trajectory_sce()
  out <- RunTrajectory(sce, reduction = "fdl", cluster_col = "seurat_clusters")
  traj <- S4Vectors::metadata(out)$SCProkaR$trajectory

  expect_equal(
    nrow(traj$paga$nodes),
    length(unique(SummarizedExperiment::colData(out)$seurat_clusters))
  )
  expect_true(all(traj$paga$edges$weight >= 0 & traj$paga$edges$weight <= 1))
})

test_that("nonzero dose panels include the shared baseline timepoint", {
  sce <- make_trajectory_sce()
  out <- RunTrajectory(
    sce,
    reduction = "fdl",
    cluster_col = "seurat_clusters",
    time_col = "time",
    dose_col = "dose",
    trajectory_col = "trajectory_state",
    baseline_dose = 0,
    baseline_time = 0
  )
  points <- S4Vectors::metadata(out)$SCProkaR$trajectory$density$points

  dose_two <- points[points$panel_dose == "2", , drop = FALSE]
  expect_true(any(as.character(dose_two$time) == "0"))
  expect_true(all(dose_two$dose[as.character(dose_two$time) == "0"] == 0))
  expect_true(all(dose_two$is_baseline[as.character(dose_two$time) == "0"]))
})

test_that("composition proportions sum to one within each dose and time", {
  sce <- make_trajectory_sce()
  out <- RunTrajectory(
    sce,
    reduction = "fdl",
    cluster_col = "seurat_clusters",
    time_col = "time",
    dose_col = "dose",
    trajectory_col = "trajectory_state"
  )
  composition <- S4Vectors::metadata(out)$SCProkaR$trajectory$composition
  grouped <- stats::aggregate(prop ~ dose + time, data = composition, FUN = sum)

  expect_true(all(abs(grouped$prop - 1) < 1e-8))
})

test_that("trajectory plotting helpers return ggplot objects", {
  sce <- make_trajectory_sce()
  out <- RunTrajectory(
    sce,
    reduction = "fdl",
    cluster_col = "seurat_clusters",
    time_col = "time",
    dose_col = "dose",
    trajectory_col = "trajectory_state"
  )

  expect_s3_class(PlotTrajectoryPAGA(out), "ggplot")
  expect_s3_class(PlotTrajectoryDensity(out), "ggplot")
  expect_s3_class(PlotTrajectoryComposition(out), "ggplot")
})
