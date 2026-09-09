## Column names that are referenced inside ggplot2::aes() and other
## data-masking contexts. Declaring them here keeps `R CMD check` from
## reporting them as undefined global variables. Every name below is a
## column of a data frame that is constructed inside the function that
## uses it; none of them is a package-level object.
utils::globalVariables(c(
    ## benchmarking.R
    "method", "metric", "score", "category_score",
    "batch_removal", "bio_conservation",
    ## integration.R
    "dim1", "dim2", ".colour",
    ## de-edger.R
    "gene", "logFC", "logCPM", "neg_log10_fdr", "de_class",
    "time", "condition", "status",
    ## tata-branch.R / tata-plot.R
    "Dim1", "Dim2", "branch", "cluster", "colour_value", "entropy",
    "direction_group", "direction_label", "edge_status",
    "fitted_expression", "gene_branch", "median", "median_time",
    "path_id", "plot_size", "plot_width", "predicted", "probability",
    "pseudotime", "q25", "q75", "signed_time_score", "tata_weight",
    "timepoint", "topology_weight", "truth", "value", "weight",
    "x", "xend", "y", "yend"
))
