toy_counts <- function() {
  matrix(
    c(
      5, 0, 1, 0, 8, 1,
      0, 3, 0, 2, 1, 0,
      9, 8, 7, 0, 0, 0,
      1, 1, 1, 9, 8, 7
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("rrsA", "geneB", "rplC", "geneD"),
      paste0("cell", 1:6)
    )
  )
}

toy_metadata <- function() {
  data.frame(
    sample_id = rep(c("sample1", "sample2"), each = 3),
    batch = rep(c("batch1", "batch2"), each = 3),
    condition = rep(c("ctrl", "treat"), each = 3),
    timepoint = c(0, 1, 2, 0, 1, 2),
    cluster = c("c1", "c1", "c2", "c2", "c3", "c3"),
    cell_type = c("A", "A", "B", "B", "C", "C"),
    row.names = paste0("cell", 1:6),
    stringsAsFactors = FALSE
  )
}

toy_feature_metadata <- function() {
  data.frame(
    feature_class = c("rrna", "gene", "ribosomal_protein", "gene"),
    row.names = c("rrsA", "geneB", "rplC", "geneD"),
    stringsAsFactors = FALSE
  )
}

make_toy_sce <- function() {
  CreateBacObject(
    toy_counts(),
    cell_metadata = toy_metadata(),
    feature_metadata = toy_feature_metadata(),
    sample_col = "sample_id",
    batch_col = "batch",
    condition_col = "condition",
    time_col = "timepoint"
  )
}

assign_toy_pseudoreplicates <- function(
    sce,
    group_cols = c("condition", "timepoint"),
    n_reps = 3,
    seed = 1) {
  set.seed(seed)

  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  strata <- interaction(meta[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  rep_id <- integer(ncol(sce))

  for (level_name in levels(strata)) {
    idx <- which(strata == level_name)
    rep_id[idx] <- sample(rep(seq_len(n_reps), length.out = length(idx)))
  }

  sce$pb_replicate <- factor(rep_id, levels = seq_len(n_reps))
  sce$pb_sample_id <- paste(
    "pb",
    as.character(sce$condition),
    paste0("t", sce$timepoint),
    paste0("r", sce$pb_replicate),
    sep = "_"
  )

  sce
}

make_toy_multidrug_sce <- function(
    n_cells = 900,
    n_features = 70,
    seed = 1) {
  simulate_tata_multidrug_sce(
    n_cells = n_cells,
    n_features = n_features,
    n_pcs = 10,
    seed = seed
  )
}

make_toy_pseudobulk <- function(
    n_cells = 900,
    n_features = 70,
    n_reps = 3,
    seed = 1) {
  sce <- make_toy_multidrug_sce(
    n_cells = n_cells,
    n_features = n_features,
    seed = seed
  )
  sce <- assign_toy_pseudoreplicates(sce, n_reps = n_reps, seed = seed)

  pb <- aggregate_pseudobulk(
    sce,
    sample_cols = c("pb_sample_id", "condition", "timepoint", "pb_replicate"),
    aggregation = "sum",
    min_cells = 10
  )

  normalize_pseudobulk(pb, method = "TMM")
}
