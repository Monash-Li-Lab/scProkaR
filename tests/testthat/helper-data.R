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
