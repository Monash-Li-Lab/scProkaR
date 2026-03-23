#' Simulate a divergent TATA dataset
#'
#' Simulate a `SingleCellExperiment` containing single-bacterium observations on
#' a divergent two-way trajectory. Time `0` is centred at the origin of the
#' latent trajectory, and later timepoints move away from that central root into
#' two branches. This makes the simulated dataset useful for demonstrating TATA
#' on a divergent trajectory where the earliest state appears in the middle of a
#' UMAP.
#'
#' The returned object contains:
#' - `counts` and `logcounts`
#' - `colData` with `cell_id`, `timepoint`, `condition`, `simulated_state`, and
#'   `simulated_branch`
#' - `reducedDims(sce)$PCA`
#' - `reducedDims(sce)$UMAP`
#' - `reducedDims(sce)$truth`
#'
#' @param n_cells Number of cells to simulate.
#' @param n_features Number of features to simulate.
#' @param n_pcs Number of principal components to store.
#' @param seed Random seed.
#'
#' @return A `SingleCellExperiment`.
#' @export
#'
#' @examples
#' sce <- simulate_tata_sce(n_cells = 500, n_features = 60, seed = 1)
#' dim(sce)
simulate_tata_sce <- function(
    n_cells = 5000,
    n_features = 100,
    n_pcs = 20,
    seed = 1) {
  simulate_tata_divergent_sce(
    n_cells = n_cells,
    n_features = n_features,
    n_pcs = n_pcs,
    seed = seed
  )
}


#' Simulate an extreme spiral-plus-divergent TATA dataset
#'
#' Simulate a stress-test `SingleCellExperiment` whose `UMAP` embedding follows
#' a spiral-like manifold that transitions into two divergent tails. This is
#' useful for visually stress-testing TATA in an intentionally extreme geometry.
#'
#' The embedding is constructed directly in UMAP-like space and paired with
#' smooth expression programs so that neighboring cells along the path remain
#' transcriptionally similar.
#'
#' @inheritParams simulate_tata_sce
#'
#' @return A `SingleCellExperiment`.
#' @export
#'
#' @examples
#' sce <- simulate_tata_spiral_sce(n_cells = 1000, n_features = 80, seed = 1)
#' SingleCellExperiment::reducedDimNames(sce)
simulate_tata_spiral_sce <- function(
    n_cells = 5000,
    n_features = 100,
    n_pcs = 20,
    seed = 1) {
  set.seed(seed)

  if (n_features < 60) {
    stop("Please simulate at least 60 features for the spiral stress test.", call. = FALSE)
  }

  branch_label <- sample(
    x = c("spiral", "tail_up", "tail_down"),
    size = n_cells,
    replace = TRUE,
    prob = c(0.72, 0.14, 0.14)
  )

  progress <- numeric(n_cells)
  progress[branch_label == "spiral"] <- stats::rbeta(sum(branch_label == "spiral"), shape1 = 1.3, shape2 = 1.5) * 0.82
  progress[branch_label == "tail_up"] <- 0.82 + stats::rbeta(sum(branch_label == "tail_up"), shape1 = 1.2, shape2 = 1.1) * 0.18
  progress[branch_label == "tail_down"] <- 0.82 + stats::rbeta(sum(branch_label == "tail_down"), shape1 = 1.2, shape2 = 1.1) * 0.18

  angle <- 1.35 * pi + 3.7 * pi * pmin(progress, 0.82) / 0.82
  radius <- 1.2 + 9.2 * pmin(progress, 0.82) / 0.82

  umap1 <- radius * cos(angle)
  umap2 <- radius * sin(angle)

  end_angle <- 1.35 * pi + 3.7 * pi
  end_radius <- 1.2 + 9.2
  end_x <- end_radius * cos(end_angle)
  end_y <- end_radius * sin(end_angle)

  tail_progress <- pmax((progress - 0.82) / 0.18, 0)
  up_idx <- branch_label == "tail_up"
  down_idx <- branch_label == "tail_down"

  if (any(up_idx)) {
    umap1[up_idx] <- end_x + 5.6 * tail_progress[up_idx] + stats::rnorm(sum(up_idx), sd = 0.40)
    umap2[up_idx] <- end_y + 0.8 + 7.5 * tail_progress[up_idx] + stats::rnorm(sum(up_idx), sd = 0.42)
  }

  if (any(down_idx)) {
    umap1[down_idx] <- end_x + 5.2 * tail_progress[down_idx] + stats::rnorm(sum(down_idx), sd = 0.42)
    umap2[down_idx] <- end_y - 0.8 - 7.2 * tail_progress[down_idx] + stats::rnorm(sum(down_idx), sd = 0.42)
  }

  spiral_noise <- 0.25 + 0.15 * (1 - pmin(progress, 0.82) / 0.82)
  spiral_idx <- branch_label == "spiral"
  umap1[spiral_idx] <- umap1[spiral_idx] + stats::rnorm(sum(spiral_idx), sd = spiral_noise[spiral_idx])
  umap2[spiral_idx] <- umap2[spiral_idx] + stats::rnorm(sum(spiral_idx), sd = spiral_noise[spiral_idx])

  timepoint <- cut(
    progress,
    breaks = c(-Inf, 0.16, 0.34, 0.54, 0.74, Inf),
    labels = c(0, 1, 2, 3, 4),
    right = TRUE
  )
  timepoint <- as.integer(as.character(timepoint))

  condition <- ifelse(
    branch_label == "tail_down",
    ifelse(stats::runif(n_cells) < 0.85, "drug", "control"),
    ifelse(branch_label == "tail_up", ifelse(stats::runif(n_cells) < 0.80, "control", "drug"), ifelse(stats::runif(n_cells) < 0.55, "control", "drug"))
  )

  state_band <- cut(
    progress,
    breaks = c(-Inf, 0.12, 0.28, 0.46, 0.66, 0.82, 0.90, Inf),
    labels = c("core", "inner", "mid", "outer", "pre_tail", "tail_early", "tail_late"),
    right = TRUE
  )

  simulated_state <- ifelse(
    branch_label == "spiral",
    paste0("spiral_", state_band),
    ifelse(branch_label == "tail_up", paste0("tail_up_", state_band), paste0("tail_down_", state_band))
  )

  cell_id <- paste0("bacterium_", seq_len(n_cells))

  meta <- data.frame(
    cell_id = cell_id,
    timepoint = timepoint,
    condition = factor(condition, levels = c("control", "drug")),
    simulated_state = factor(simulated_state),
    simulated_branch = factor(branch_label, levels = c("spiral", "tail_up", "tail_down")),
    progress = progress,
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$cell_id

  progress_signal <- progress
  branch_up_signal <- as.numeric(branch_label == "tail_up") * tail_progress
  branch_down_signal <- as.numeric(branch_label == "tail_down") * tail_progress
  radius_signal <- pmin(radius / max(radius), 1)
  angle_signal <- sin(angle)
  angle_cos_signal <- cos(angle)
  root_signal <- exp(-((progress - 0.06) ^ 2) / (2 * 0.06 ^ 2))
  inner_signal <- exp(-((progress - 0.20) ^ 2) / (2 * 0.08 ^ 2))
  mid_signal <- exp(-((progress - 0.42) ^ 2) / (2 * 0.09 ^ 2))
  outer_signal <- exp(-((progress - 0.66) ^ 2) / (2 * 0.08 ^ 2))
  tail_signal <- pmax((progress - 0.78) / 0.22, 0)
  condition_signal <- 0.20 * as.numeric(condition == "drug") + 0.22 * branch_down_signal

  module_names <- c(
    "root",
    "inner",
    "mid",
    "outer",
    "progress",
    "angle_sin",
    "angle_cos",
    "tail_up",
    "tail_down",
    "condition",
    "noise"
  )

  module_sizes <- round(c(0.10, 0.10, 0.12, 0.10, 0.16, 0.08, 0.08, 0.08, 0.08, 0.06, 0.04) * n_features)
  module_sizes[1] <- module_sizes[1] + (n_features - sum(module_sizes))
  feature_module <- rep(module_names, times = module_sizes)

  baseline <- stats::rnorm(n_features, mean = 0.78, sd = 0.12)
  root_load <- rep(0, n_features)
  inner_load <- rep(0, n_features)
  mid_load <- rep(0, n_features)
  outer_load <- rep(0, n_features)
  progress_load <- rep(0, n_features)
  angle_sin_load <- rep(0, n_features)
  angle_cos_load <- rep(0, n_features)
  tail_up_load <- rep(0, n_features)
  tail_down_load <- rep(0, n_features)
  condition_load <- rep(0, n_features)

  root_load[feature_module == "root"] <- stats::runif(sum(feature_module == "root"), 0.70, 0.95)
  inner_load[feature_module == "inner"] <- stats::runif(sum(feature_module == "inner"), 0.60, 0.85)
  mid_load[feature_module == "mid"] <- stats::runif(sum(feature_module == "mid"), 0.60, 0.85)
  outer_load[feature_module == "outer"] <- stats::runif(sum(feature_module == "outer"), 0.60, 0.85)
  progress_load[feature_module == "progress"] <- stats::runif(sum(feature_module == "progress"), 0.60, 0.90)
  angle_sin_load[feature_module == "angle_sin"] <- stats::runif(sum(feature_module == "angle_sin"), 0.35, 0.55)
  angle_cos_load[feature_module == "angle_cos"] <- stats::runif(sum(feature_module == "angle_cos"), 0.35, 0.55)
  tail_up_load[feature_module == "tail_up"] <- stats::runif(sum(feature_module == "tail_up"), 0.55, 0.85)
  tail_down_load[feature_module == "tail_down"] <- stats::runif(sum(feature_module == "tail_down"), 0.55, 0.85)
  condition_load[feature_module == "condition"] <- stats::runif(sum(feature_module == "condition"), 0.18, 0.32)

  eta <- matrix(rep(baseline, each = n_cells), nrow = n_cells) +
    outer(root_signal, root_load) +
    outer(inner_signal, inner_load) +
    outer(mid_signal, mid_load) +
    outer(outer_signal + tail_signal, outer_load) +
    outer(progress_signal + radius_signal, progress_load) +
    outer(angle_signal, angle_sin_load) +
    outer(angle_cos_signal, angle_cos_load) +
    outer(branch_up_signal, tail_up_load) +
    outer(branch_down_signal, tail_down_load) +
    outer(condition_signal, condition_load) +
    matrix(stats::rnorm(n_cells * n_features, mean = 0, sd = 0.07), nrow = n_cells)

  eta <- pmax(pmin(eta, 4.2), -1.5)
  mu <- exp(eta)
  mu <- t(mu)

  size_factor <- stats::rgamma(n_cells, shape = 10, rate = 10)
  mu <- sweep(mu, 2, size_factor * 6.5, `*`)

  counts <- matrix(
    stats::rnbinom(n = length(mu), mu = as.vector(mu), size = 10),
    nrow = n_features,
    ncol = n_cells
  )

  gene_ids <- paste0("gene_", seq_len(n_features))
  colnames(counts) <- meta$cell_id
  rownames(counts) <- gene_ids

  library_size <- colSums(counts)
  norm_factor <- library_size / stats::median(library_size)
  logcounts <- log1p(sweep(counts, 2, norm_factor, `/`))
  rownames(logcounts) <- gene_ids
  colnames(logcounts) <- meta$cell_id

  latent_for_pca <- cbind(
    progress = progress_signal,
    radius = radius_signal,
    angle_sin = angle_signal,
    angle_cos = angle_cos_signal,
    tail_up = branch_up_signal,
    tail_down = branch_down_signal,
    condition = condition_signal,
    umap1 = scale(umap1)[, 1],
    umap2 = scale(umap2)[, 1]
  )
  pca <- stats::prcomp(latent_for_pca, center = TRUE, scale. = TRUE)
  n_keep <- min(n_pcs, ncol(pca$x))
  pcs <- pca$x[, seq_len(n_keep), drop = FALSE]
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  rownames(pcs) <- meta$cell_id

  truth <- cbind(latent_progress = progress, branch_code = c(spiral = 0, tail_up = 1, tail_down = -1)[branch_label])
  rownames(truth) <- meta$cell_id

  umap <- cbind(UMAP1 = umap1, UMAP2 = umap2)
  rownames(umap) <- meta$cell_id

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = logcounts
    ),
    colData = S4Vectors::DataFrame(meta)
  )

  SingleCellExperiment::reducedDim(sce, "PCA") <- pcs
  SingleCellExperiment::reducedDim(sce, "UMAP") <- umap
  SingleCellExperiment::reducedDim(sce, "truth") <- truth
  SummarizedExperiment::rowData(sce)$feature_module <- feature_module

  S4Vectors::metadata(sce)$simulation <- list(
    seed = seed,
    trajectory = "spiral_divergent",
    feature_module = feature_module
  )

  sce
}


#' Simulate a multi-drug branching TATA dataset
#'
#' Simulate a `SingleCellExperiment` with one shared origin trajectory that can
#' diverge into three drug-associated fates. The three drug-associated branches
#' emerge at different real treatment times, making the dataset useful for
#' checking whether TATA uses time information rather than feature similarity
#' alone when orienting cluster-level transitions.
#'
#' The simulated design uses treatment times in hours:
#' - `drug_A` branch emerges early, around 24 hours
#' - `drug_C` branch emerges around 48 hours
#' - `drug_B` branch emerges late, around 96 hours
#'
#' The geometry is intentionally center-rooted rather than tip-rooted:
#' `drug_A` and `drug_B` extend in opposite directions from a shared middle
#' state, and `drug_C` diverges along a third arm. This is meant to stress-test
#' whether time labels can recover direction when the root is not visually
#' obvious from topology alone.
#'
#' @inheritParams simulate_tata_sce
#'
#' @return A `SingleCellExperiment`.
#' @export
#'
#' @examples
#' sce <- simulate_tata_multidrug_sce(n_cells = 1000, n_features = 80, seed = 1)
#' table(sce$condition, sce$timepoint)
simulate_tata_multidrug_sce <- function(
    n_cells = 5000,
    n_features = 100,
    n_pcs = 20,
    seed = 1) {
  set.seed(seed)

  if (n_features < 60) {
    stop("Please simulate at least 60 features for the multi-drug branching example.", call. = FALSE)
  }

  condition_levels <- c("drug_A", "drug_B", "drug_C")
  time_levels <- c(0, 24, 48, 72, 96, 120)
  branch_onset <- c(drug_A = 24, drug_B = 96, drug_C = 48)

  design <- expand.grid(
    condition = condition_levels,
    timepoint = time_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  design$weight <- c(
    420, 520, 480, 430, 380, 320,
    420, 520, 500, 470, 430, 360,
    420, 520, 500, 430, 360, 320
  )
  design$n <- round(n_cells * design$weight / sum(design$weight))
  delta <- n_cells - sum(design$n)
  if (delta != 0L) {
    adjust_idx <- seq_len(abs(delta))
    design$n[adjust_idx] <- design$n[adjust_idx] + sign(delta)
  }

  meta <- design[rep(seq_len(nrow(design)), design$n), c("condition", "timepoint"), drop = FALSE]
  rownames(meta) <- NULL
  meta <- meta[sample(seq_len(nrow(meta))), , drop = FALSE]

  cell_id <- paste0("bacterium_", seq_len(nrow(meta)))
  meta$cell_id <- cell_id

  time_scaled <- meta$timepoint / max(time_levels)
  onset_scaled <- branch_onset[meta$condition] / max(time_levels)

  # Continuous progression within each labeled timepoint creates many
  # transitional cells instead of isolated condition-time islands.
  progress <- pmin(
    1,
    pmax(
      0,
      time_scaled + stats::rnorm(nrow(meta), mean = 0, sd = 0.07)
    )
  )

  branch_activation <- pmax(
    0,
    progress - onset_scaled + stats::rnorm(nrow(meta), mean = 0, sd = 0.020)
  )
  branch_activation <- pmin(0.48, branch_activation) / 0.48
  branch_activation <- pmin(1, pmax(0, branch_activation))

  branch_direction <- rbind(
    drug_A = c(-1.18, 0.12),
    drug_B = c(1.18, -0.10),
    drug_C = c(0.05, 1.22)
  )

  dir_x <- branch_direction[meta$condition, 1]
  dir_y <- branch_direction[meta$condition, 2]

  central_radius <- 0.25 + 1.05 * pmin(progress, 0.35) / 0.35
  central_angle <- 1.6 * pi * progress + stats::rnorm(nrow(meta), mean = 0, sd = 0.08)
  trunk_x <- central_radius * cos(central_angle) * (1 - branch_activation)
  trunk_y <- -0.8 + 1.8 * progress + 0.7 * sin(1.2 * pi * progress) * (1 - 0.4 * branch_activation)

  branch_length <- 7.0 * branch_activation ^ 1.10
  branch_bend <- 1.2 * branch_activation * (progress - onset_scaled)
  bend_sign <- c(drug_A = 1, drug_B = -1, drug_C = 1)[meta$condition]

  umap1 <- trunk_x +
    branch_length * dir_x +
    0.9 * branch_bend * ifelse(meta$condition == "drug_C", 0.35, 1) +
    stats::rnorm(nrow(meta), sd = 0.18 + 0.16 * (1 - branch_activation))

  umap2 <- trunk_y +
    branch_length * dir_y +
    0.8 * bend_sign * branch_bend * ifelse(meta$condition == "drug_C", 0.15, 0.55) +
    stats::rnorm(nrow(meta), sd = 0.18 + 0.16 * (1 - branch_activation))

  fate_label <- ifelse(
    branch_activation < 0.10,
    "origin",
    ifelse(branch_activation < 0.55, "transition", "fate")
  )
  simulated_state <- paste(meta$condition, fate_label, paste0("h", meta$timepoint), sep = "_")

  meta$condition <- factor(meta$condition, levels = condition_levels)
  meta$timepoint <- as.numeric(meta$timepoint)
  meta$simulated_state <- factor(simulated_state)
  meta$simulated_branch <- factor(
    ifelse(branch_activation < 0.10, "origin", as.character(meta$condition)),
    levels = c("origin", condition_levels)
  )
  meta$branch_activation <- branch_activation
  meta$progress <- progress
  meta <- meta[, c("cell_id", "timepoint", "condition", "simulated_state", "simulated_branch", "progress", "branch_activation")]
  rownames(meta) <- meta$cell_id

  root_signal <- exp(-((progress - 0.08) ^ 2) / (2 * 0.06 ^ 2))
  early_signal <- exp(-((progress - 0.24) ^ 2) / (2 * 0.08 ^ 2))
  mid_signal <- exp(-((progress - 0.50) ^ 2) / (2 * 0.11 ^ 2))
  late_signal <- exp(-((progress - 0.78) ^ 2) / (2 * 0.12 ^ 2))
  center_signal <- exp(-((progress - 0.14) ^ 2) / (2 * 0.10 ^ 2)) * (1 - branch_activation)
  branch_A_signal <- branch_activation * as.numeric(meta$condition == "drug_A")
  branch_B_signal <- branch_activation * as.numeric(meta$condition == "drug_B")
  branch_C_signal <- branch_activation * as.numeric(meta$condition == "drug_C")
  transition_signal <- pmax(0, 1 - abs(branch_activation - 0.28) / 0.28)
  # Condition-specific offsets should only appear once a branch is activated;
  # otherwise the early shared origin would split by treatment too soon.
  condition_signal <- branch_activation * (
    0.12 * as.numeric(meta$condition == "drug_B") +
      0.08 * as.numeric(meta$condition == "drug_C")
  )
  umap1_signal <- scale(umap1)[, 1]
  umap2_signal <- scale(umap2)[, 1]

  module_names <- c(
    "root",
    "early",
    "mid",
    "late",
    "center",
    "transition",
    "branch_A",
    "branch_B",
    "branch_C",
    "condition",
    "umap1",
    "umap2",
    "noise"
  )

  module_sizes <- round(c(0.09, 0.09, 0.10, 0.10, 0.08, 0.10, 0.09, 0.09, 0.09, 0.07, 0.06, 0.07, 0.03) * n_features)
  module_sizes[1] <- module_sizes[1] + (n_features - sum(module_sizes))
  feature_module <- rep(module_names, times = module_sizes)

  baseline <- stats::rnorm(n_features, mean = 0.82, sd = 0.13)
  root_load <- rep(0, n_features)
  early_load <- rep(0, n_features)
  mid_load <- rep(0, n_features)
  late_load <- rep(0, n_features)
  center_load <- rep(0, n_features)
  transition_load <- rep(0, n_features)
  branch_A_load <- rep(0, n_features)
  branch_B_load <- rep(0, n_features)
  branch_C_load <- rep(0, n_features)
  condition_load <- rep(0, n_features)
  umap1_load <- rep(0, n_features)
  umap2_load <- rep(0, n_features)

  root_load[feature_module == "root"] <- stats::runif(sum(feature_module == "root"), 0.75, 1.00)
  early_load[feature_module == "early"] <- stats::runif(sum(feature_module == "early"), 0.60, 0.88)
  mid_load[feature_module == "mid"] <- stats::runif(sum(feature_module == "mid"), 0.60, 0.88)
  late_load[feature_module == "late"] <- stats::runif(sum(feature_module == "late"), 0.60, 0.88)
  center_load[feature_module == "center"] <- stats::runif(sum(feature_module == "center"), 0.40, 0.65)
  transition_load[feature_module == "transition"] <- stats::runif(sum(feature_module == "transition"), 0.45, 0.70)
  branch_A_load[feature_module == "branch_A"] <- stats::runif(sum(feature_module == "branch_A"), 0.45, 0.72)
  branch_B_load[feature_module == "branch_B"] <- stats::runif(sum(feature_module == "branch_B"), 0.45, 0.72)
  branch_C_load[feature_module == "branch_C"] <- stats::runif(sum(feature_module == "branch_C"), 0.45, 0.72)
  condition_load[feature_module == "condition"] <- stats::runif(sum(feature_module == "condition"), 0.16, 0.30)
  umap1_load[feature_module == "umap1"] <- stats::runif(sum(feature_module == "umap1"), 0.18, 0.35)
  umap2_load[feature_module == "umap2"] <- stats::runif(sum(feature_module == "umap2"), 0.18, 0.35)

  eta <- matrix(rep(baseline, each = nrow(meta)), nrow = nrow(meta)) +
    outer(root_signal, root_load) +
    outer(early_signal, early_load) +
    outer(mid_signal, mid_load) +
    outer(late_signal, late_load) +
    outer(center_signal, center_load) +
    outer(transition_signal, transition_load) +
    outer(branch_A_signal, branch_A_load) +
    outer(branch_B_signal, branch_B_load) +
    outer(branch_C_signal, branch_C_load) +
    outer(condition_signal, condition_load) +
    outer(umap1_signal, umap1_load) +
    outer(umap2_signal, umap2_load) +
    matrix(stats::rnorm(nrow(meta) * n_features, mean = 0, sd = 0.08), nrow = nrow(meta))

  eta <- pmax(pmin(eta, 4.2), -1.5)
  mu <- exp(eta)
  mu <- t(mu)

  size_factor <- stats::rgamma(nrow(meta), shape = 10, rate = 10)
  mu <- sweep(mu, 2, size_factor * 6.8, `*`)

  counts <- matrix(
    stats::rnbinom(n = length(mu), mu = as.vector(mu), size = 10),
    nrow = n_features,
    ncol = nrow(meta)
  )

  gene_ids <- paste0("gene_", seq_len(n_features))
  colnames(counts) <- meta$cell_id
  rownames(counts) <- gene_ids

  library_size <- colSums(counts)
  norm_factor <- library_size / stats::median(library_size)
  logcounts <- log1p(sweep(counts, 2, norm_factor, `/`))
  rownames(logcounts) <- gene_ids
  colnames(logcounts) <- meta$cell_id

  latent_for_pca <- cbind(
    progress = progress,
    branch_A = branch_A_signal,
    branch_B = branch_B_signal,
    branch_C = branch_C_signal,
    transition = transition_signal,
    condition = condition_signal,
    umap1 = umap1_signal,
    umap2 = umap2_signal
  )
  pca <- stats::prcomp(latent_for_pca, center = TRUE, scale. = TRUE)
  n_keep <- min(n_pcs, ncol(pca$x))
  pcs <- pca$x[, seq_len(n_keep), drop = FALSE]
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  rownames(pcs) <- meta$cell_id

  truth <- cbind(
    progress = progress,
    branch_activation = branch_activation,
    branch_code = c(origin = 0, drug_A = 1, drug_B = 2, drug_C = 3)[as.character(meta$simulated_branch)]
  )
  rownames(truth) <- meta$cell_id

  umap <- cbind(UMAP1 = umap1, UMAP2 = umap2)
  rownames(umap) <- meta$cell_id

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = logcounts
    ),
    colData = S4Vectors::DataFrame(meta)
  )

  SingleCellExperiment::reducedDim(sce, "PCA") <- pcs
  SingleCellExperiment::reducedDim(sce, "UMAP") <- umap
  SingleCellExperiment::reducedDim(sce, "truth") <- truth
  SummarizedExperiment::rowData(sce)$feature_module <- feature_module

  S4Vectors::metadata(sce)$simulation <- list(
    seed = seed,
    trajectory = "multidrug_branching",
    branch_onset_hours = branch_onset,
    feature_module = feature_module
  )

  sce
}


#' @rdname simulate_tata_sce
#' @export
simulate_tata_divergent_sce <- function(
    n_cells = 5000,
    n_features = 100,
    n_pcs = 20,
    seed = 1) {
  set.seed(seed)

  if (n_cells < 1000) {
    warning(
      "The divergent simulation is designed for larger datasets; ",
      "consider at least 1000 cells for clearer branch structure.",
      call. = FALSE
    )
  }

  if (n_features < 50) {
    stop("Please simulate at least 50 features for this example.", call. = FALSE)
  }

  allocation <- data.frame(
    timepoint = c(0, 1, 1, 2, 2, 3, 3, 4, 4),
    simulated_branch = c("root", "adaptive", "collapse", "adaptive", "collapse", "adaptive", "collapse", "adaptive", "collapse"),
    weight = c(900, 700, 700, 600, 600, 500, 500, 250, 250),
    stringsAsFactors = FALSE
  )

  allocation$n <- round(n_cells * allocation$weight / sum(allocation$weight))
  delta <- n_cells - sum(allocation$n)
  if (delta != 0L) {
    adjust_idx <- seq_len(abs(delta))
    allocation$n[adjust_idx] <- allocation$n[adjust_idx] + sign(delta)
  }

  meta <- allocation[rep(seq_len(nrow(allocation)), allocation$n), c("timepoint", "simulated_branch"), drop = FALSE]
  meta <- meta[sample(seq_len(nrow(meta))), , drop = FALSE]
  rownames(meta) <- NULL

  control_probability <- ifelse(
    meta$simulated_branch == "adaptive",
    0.75,
    ifelse(meta$simulated_branch == "collapse", 0.25, 0.50)
  )

  meta$condition <- ifelse(stats::runif(nrow(meta)) < control_probability, "control", "drug")
  meta$simulated_state <- ifelse(
    meta$simulated_branch == "root",
    "root",
    paste(meta$simulated_branch, paste0("t", meta$timepoint), sep = "_")
  )
  meta$cell_id <- paste0("bacterium_", seq_len(nrow(meta)))

  state_levels <- c(
    "root",
    "adaptive_t1", "adaptive_t2", "adaptive_t3", "adaptive_t4",
    "collapse_t1", "collapse_t2", "collapse_t3", "collapse_t4"
  )

  meta$condition <- factor(meta$condition, levels = c("control", "drug"))
  meta$simulated_branch <- factor(meta$simulated_branch, levels = c("root", "adaptive", "collapse"))
  meta$simulated_state <- factor(meta$simulated_state, levels = state_levels)
  meta <- meta[, c("cell_id", "timepoint", "condition", "simulated_state", "simulated_branch")]
  rownames(meta) <- meta$cell_id

  branch_sign <- ifelse(meta$simulated_branch == "adaptive", -1, ifelse(meta$simulated_branch == "collapse", 1, 0))
  is_drug <- as.numeric(meta$condition == "drug")

  # Each experimental timepoint occupies a broad, overlapping interval on a
  # continuous latent progression axis. This creates many intermediate cells
  # between labelled timepoints, which helps UMAP recover a trajectory rather
  # than disconnected blobs.
  tau_min <- c(`0` = 0.00, `1` = 0.08, `2` = 0.24, `3` = 0.44, `4` = 0.64)
  tau_max <- c(`0` = 0.18, `1` = 0.34, `2` = 0.58, `3` = 0.80, `4` = 1.00)
  tau_continuous <- stats::runif(
    nrow(meta),
    min = tau_min[as.character(meta$timepoint)],
    max = tau_max[as.character(meta$timepoint)]
  )
  tau_continuous <- pmin(1, pmax(0, tau_continuous + stats::rnorm(nrow(meta), mean = 0, sd = 0.025)))

  # Branch commitment turns on gradually after the root. Near the split, cells
  # from the two future branches still lie close together, generating visible
  # transition states in PCA/UMAP space.
  branch_strength <- ifelse(
    branch_sign == 0,
    0,
    plogis((tau_continuous - 0.34) / 0.09) + stats::rnorm(nrow(meta), mean = 0, sd = 0.035)
  )
  branch_strength <- pmin(1, pmax(0, branch_strength))
  branch_effect <- branch_sign * branch_strength

  latent_x <- 4.20 * branch_effect * (0.55 + 0.55 * tau_continuous) +
    stats::rnorm(nrow(meta), mean = 0, sd = 0.11 + 0.05 * (1 - branch_strength))

  latent_y <- 5.00 * tau_continuous - 1.15 +
    0.45 * tau_continuous^2 -
    0.15 * branch_strength^2 +
    stats::rnorm(nrow(meta), mean = 0, sd = 0.10)

  root_signal <- exp(-((tau_continuous - 0.05) ^ 2) / (2 * 0.08 ^ 2))
  early_signal <- exp(-((tau_continuous - 0.22) ^ 2) / (2 * 0.10 ^ 2))
  mid_signal <- exp(-((tau_continuous - 0.48) ^ 2) / (2 * 0.14 ^ 2))
  late_signal <- exp(-((tau_continuous - 0.78) ^ 2) / (2 * 0.14 ^ 2))
  branch_signal <- branch_effect
  adaptive_signal <- pmax(-branch_signal, 0)
  collapse_signal <- pmax(branch_signal, 0)
  stress_signal <- 0.25 * is_drug + 0.30 * collapse_signal + 0.10 * tau_continuous
  condition_signal <- 0.18 * is_drug + 0.08 * tau_continuous * is_drug

  module_names <- c(
    "root",
    "early",
    "mid",
    "late",
    "branch",
    "adaptive",
    "collapse",
    "stress",
    "latent_y",
    "noise"
  )

  module_sizes <- round(c(0.14, 0.14, 0.16, 0.12, 0.12, 0.08, 0.08, 0.07, 0.05, 0.04) * n_features)
  module_sizes[1] <- module_sizes[1] + (n_features - sum(module_sizes))
  feature_module <- rep(module_names, times = module_sizes)

  baseline <- stats::rnorm(n_features, mean = 0.82, sd = 0.14)
  root_load <- rep(0, n_features)
  early_load <- rep(0, n_features)
  mid_load <- rep(0, n_features)
  late_load <- rep(0, n_features)
  branch_load <- rep(0, n_features)
  adaptive_load <- rep(0, n_features)
  collapse_load <- rep(0, n_features)
  stress_load <- rep(0, n_features)
  latent_y_load <- rep(0, n_features)

  root_load[feature_module == "root"] <- stats::runif(sum(feature_module == "root"), 0.75, 1.05)
  early_load[feature_module == "early"] <- stats::runif(sum(feature_module == "early"), 0.65, 0.95)
  mid_load[feature_module == "mid"] <- stats::runif(sum(feature_module == "mid"), 0.65, 0.95)
  late_load[feature_module == "late"] <- stats::runif(sum(feature_module == "late"), 0.65, 0.95)
  branch_load[feature_module == "branch"] <- rep(c(-0.50, 0.50), length.out = sum(feature_module == "branch"))
  adaptive_load[feature_module == "adaptive"] <- stats::runif(sum(feature_module == "adaptive"), 0.30, 0.55)
  collapse_load[feature_module == "collapse"] <- stats::runif(sum(feature_module == "collapse"), 0.30, 0.55)
  stress_load[feature_module == "stress"] <- stats::runif(sum(feature_module == "stress"), 0.18, 0.35)
  latent_y_load[feature_module == "latent_y"] <- stats::runif(sum(feature_module == "latent_y"), 0.22, 0.40)

  eta <- matrix(rep(baseline, each = nrow(meta)), nrow = nrow(meta)) +
    outer(root_signal, root_load) +
    outer(early_signal, early_load) +
    outer(mid_signal, mid_load) +
    outer(late_signal, late_load) +
    outer(branch_signal, branch_load) +
    outer(adaptive_signal, adaptive_load) +
    outer(collapse_signal, collapse_load) +
    outer(stress_signal + condition_signal, stress_load) +
    outer(latent_y, latent_y_load) +
    matrix(stats::rnorm(nrow(meta) * n_features, mean = 0, sd = 0.08), nrow = nrow(meta))

  eta <- pmax(pmin(eta, 4.5), -1.5)
  mu <- exp(eta)
  mu <- t(mu)

  size_factor <- stats::rgamma(nrow(meta), shape = 10, rate = 10)
  mu <- sweep(mu, 2, size_factor * 7, `*`)

  counts <- matrix(
    stats::rnbinom(n = length(mu), mu = as.vector(mu), size = 8),
    nrow = n_features,
    ncol = nrow(meta)
  )

  gene_ids <- paste0("gene_", seq_len(n_features))
  colnames(counts) <- meta$cell_id
  rownames(counts) <- gene_ids

  library_size <- colSums(counts)
  norm_factor <- library_size / stats::median(library_size)
  logcounts <- log1p(sweep(counts, 2, norm_factor, `/`))
  rownames(logcounts) <- gene_ids
  colnames(logcounts) <- meta$cell_id

  pca <- stats::prcomp(t(logcounts), center = TRUE, scale. = TRUE)
  n_keep <- min(n_pcs, ncol(pca$x))
  pcs <- pca$x[, seq_len(n_keep), drop = FALSE]
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  rownames(pcs) <- meta$cell_id

  truth <- cbind(latent_x = latent_x, latent_y = latent_y)
  rownames(truth) <- meta$cell_id

  umap <- cbind(
    UMAP1 = latent_x + 0.15 * sin(latent_y),
    UMAP2 = 0.85 * latent_y + 0.10 * sin(latent_x / 2)
  )
  rownames(umap) <- meta$cell_id

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = logcounts
    ),
    colData = S4Vectors::DataFrame(meta)
  )

  SingleCellExperiment::reducedDim(sce, "PCA") <- pcs
  SingleCellExperiment::reducedDim(sce, "UMAP") <- umap
  SingleCellExperiment::reducedDim(sce, "truth") <- truth
  SummarizedExperiment::rowData(sce)$feature_module <- feature_module

  S4Vectors::metadata(sce)$simulation <- list(
    allocation = allocation,
    feature_module = feature_module,
    seed = seed,
    trajectory = "divergent"
  )

  sce
}
