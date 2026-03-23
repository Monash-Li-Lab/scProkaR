#' Plot the cluster-level TATA graph
#'
#' Plot the directed cluster-level TATA graph. The graph can be drawn in its
#' own abstract layout or positioned on a chosen reduced dimension so that the
#' graph coordinates match the embedding shown to the user.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param dimred Optional reduced dimension used to position cluster nodes. When
#'   supplied together with `layout_mode = "embedding"`, cluster centroids are
#'   computed from this embedding.
#' @param layout_mode Layout mode for the graph. Choices are:
#' - `"graph"`: use the stored abstract graph coordinates, or a force-directed
#'   layout if they are unavailable.
#' - `"embedding"`: place cluster nodes at their centroids in the selected
#'   reduced dimension.
#' @param cluster_col Optional cluster column used when `layout_mode =
#'   "embedding"`. By default this is taken from `tata_result$parameters`.
#' @param show_cells Logical indicating whether to draw cells in the background
#'   when `layout_mode = "embedding"`.
#' @param colour_by Optional `colData` column used to color background cells. If
#'   `NULL` and `show_cells = TRUE`, the cluster column is used when available.
#' @param node_colour_by Node coloring scheme. Choices are:
#' - `"median_time"`: color nodes by cluster median time.
#' - `"cluster"`: color nodes by cluster identity, which is useful for
#'   embedding-aligned graph plots.
#' @param main Plot title.
#' @param vertex_label_size Label size.
#' @param edge_width_scale Relative edge width scale.
#' @param point_size Background point size when cells are drawn.
#' @param point_alpha Background point alpha when cells are drawn.
#' @param arrow_size Arrow size for directed edges.
#'
#' @return A `ggplot2` object.
#' @export
plot_tata_cluster_graph <- function(
    tata_result,
    dimred = NULL,
    layout_mode = NULL,
    cluster_col = NULL,
    show_cells = FALSE,
    colour_by = NULL,
    node_colour_by = c("median_time", "cluster"),
    main = "TATA cluster graph",
    vertex_label_size = 3.5,
    edge_width_scale = 5,
    point_size = 0.4,
    point_alpha = 0.45,
    arrow_size = 0.15) {
  cluster_graph <- tata_result$cluster_graph

  if (!inherits(cluster_graph, "igraph")) {
    stop("`tata_result$cluster_graph` must be an igraph object.", call. = FALSE)
  }

  if (igraph::vcount(cluster_graph) == 0L) {
    stop("The cluster graph is empty.", call. = FALSE)
  }

  if (is.null(layout_mode)) {
    layout_mode <- if (!is.null(dimred)) "embedding" else "graph"
  }
  layout_mode <- match.arg(layout_mode, c("graph", "embedding"))
  node_colour_by <- match.arg(node_colour_by)

  graph_vertex_df <- data.frame(
    cluster = igraph::V(cluster_graph)$name,
    size = igraph::V(cluster_graph)$size,
    median_time = igraph::V(cluster_graph)$median_time,
    stringsAsFactors = FALSE
  )

  cell_df <- NULL
  if (identical(layout_mode, "embedding")) {
    sce <- tata_result$sce
    if (!methods::is(sce, "SingleCellExperiment")) {
      stop(
        "`layout_mode = \"embedding\"` requires `tata_result$sce` to be a SingleCellExperiment.",
        call. = FALSE
      )
    }

    if (is.null(dimred)) {
      if ("UMAP" %in% SingleCellExperiment::reducedDimNames(sce)) {
        dimred <- "UMAP"
      } else {
        stop("Please supply `dimred` for `layout_mode = \"embedding\"`.", call. = FALSE)
      }
    }

    if (!dimred %in% SingleCellExperiment::reducedDimNames(sce)) {
      stop("Reduced dimension `", dimred, "` is not present in `tata_result$sce`.", call. = FALSE)
    }

    if (is.null(cluster_col)) {
      cluster_col <- tata_result$parameters$cluster_col
    }

    if (is.null(cluster_col) || !cluster_col %in% colnames(SummarizedExperiment::colData(sce))) {
      stop("A valid `cluster_col` is required for `layout_mode = \"embedding\"`.", call. = FALSE)
    }

    embedding <- as.matrix(SingleCellExperiment::reducedDim(sce, dimred))
    if (ncol(embedding) < 2L) {
      stop("The requested reduced dimension must contain at least two columns.", call. = FALSE)
    }

    meta <- as.data.frame(SummarizedExperiment::colData(sce))
    meta$cluster_plot <- as.character(meta[[cluster_col]])
    meta$Dim1 <- embedding[, 1]
    meta$Dim2 <- embedding[, 2]

    node_df <- stats::aggregate(
      meta[, c("Dim1", "Dim2")],
      by = list(cluster = meta$cluster_plot),
      FUN = stats::median
    )
    names(node_df)[2:3] <- c("x", "y")
    node_df <- merge(graph_vertex_df, node_df, by = "cluster", all.x = TRUE, sort = FALSE)
    node_df <- node_df[match(graph_vertex_df$cluster, node_df$cluster), , drop = FALSE]

    if (isTRUE(show_cells)) {
      if (is.null(colour_by)) {
        colour_by <- cluster_col
      }

      cell_df <- meta[, c("Dim1", "Dim2", "cluster_plot"), drop = FALSE]
      if (!is.null(colour_by) && colour_by %in% colnames(meta)) {
        cell_df$colour_value <- meta[[colour_by]]
      } else {
        cell_df$colour_value <- "cells"
      }
    }
  } else {
    if (isTRUE(show_cells)) {
      warning("`show_cells = TRUE` is only used with `layout_mode = \"embedding\"`.", call. = FALSE)
    }

    layout_mat <- cbind(igraph::V(cluster_graph)$x, igraph::V(cluster_graph)$y)
    if (any(!is.finite(layout_mat))) {
      layout_mat <- igraph::layout_with_fr(cluster_graph)
    }

    node_df <- graph_vertex_df
    node_df$x <- layout_mat[, 1]
    node_df$y <- layout_mat[, 2]
  }

  if (any(!is.finite(node_df$x)) || any(!is.finite(node_df$y))) {
    stop("Cluster node positions could not be computed for plotting.", call. = FALSE)
  }

  edge_df <- igraph::as_data_frame(cluster_graph, what = "edges")
  if (nrow(edge_df) > 0L) {
    from_idx <- match(edge_df$from, node_df$cluster)
    to_idx <- match(edge_df$to, node_df$cluster)
    edge_df$x <- node_df$x[from_idx]
    edge_df$y <- node_df$y[from_idx]
    edge_df$xend <- node_df$x[to_idx]
    edge_df$yend <- node_df$y[to_idx]
  } else {
    edge_df <- data.frame(
      from = character(0),
      to = character(0),
      weight = numeric(0),
      ambiguous = logical(0),
      x = numeric(0),
      y = numeric(0),
      xend = numeric(0),
      yend = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  if (!"ambiguous" %in% colnames(edge_df)) {
    edge_df$ambiguous <- FALSE
  }

  if (nrow(edge_df) > 0L && max(edge_df$weight, na.rm = TRUE) > 0) {
    edge_df$plot_width <- 0.5 + edge_width_scale * edge_df$weight / max(edge_df$weight, na.rm = TRUE)
  } else {
    edge_df$plot_width <- rep(0.5, nrow(edge_df))
  }

  if (length(unique(node_df$size)) > 1L) {
    node_df$plot_size <- 4 + 6 * (node_df$size - min(node_df$size)) / diff(range(node_df$size))
  } else {
    node_df$plot_size <- rep(7, nrow(node_df))
  }

  node_df$cluster <- factor(node_df$cluster, levels = node_df$cluster)
  cluster_palette <- grDevices::hcl.colors(nrow(node_df), palette = "Dark 3")
  names(cluster_palette) <- as.character(node_df$cluster)

  p <- ggplot2::ggplot()

  if (!is.null(cell_df) && nrow(cell_df) > 0L) {
    if (!is.null(colour_by) && colour_by == cluster_col) {
      cell_df$colour_value <- factor(as.character(cell_df$colour_value), levels = levels(node_df$cluster))
      p <- p + ggplot2::geom_point(
        data = cell_df,
        mapping = ggplot2::aes(x = Dim1, y = Dim2, colour = colour_value),
        size = point_size,
        alpha = point_alpha
      ) +
        ggplot2::scale_colour_manual(values = cluster_palette, name = colour_by)
    } else if (!is.null(colour_by) && is.numeric(cell_df$colour_value)) {
      p <- p + ggplot2::geom_point(
        data = cell_df,
        mapping = ggplot2::aes(x = Dim1, y = Dim2, colour = colour_value),
        size = point_size,
        alpha = point_alpha
      ) +
        ggplot2::scale_colour_viridis_c(option = "plasma", name = colour_by)
    } else if (!is.null(colour_by) && !is.numeric(cell_df$colour_value)) {
      p <- p + ggplot2::geom_point(
        data = cell_df,
        mapping = ggplot2::aes(x = Dim1, y = Dim2, colour = colour_value),
        size = point_size,
        alpha = point_alpha
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = cell_df,
        mapping = ggplot2::aes(x = Dim1, y = Dim2),
        colour = "grey80",
        size = point_size,
        alpha = point_alpha
      )
    }
  }

  directed_edges <- edge_df[!edge_df$ambiguous, , drop = FALSE]
  ambiguous_edges <- edge_df[edge_df$ambiguous, , drop = FALSE]

  if (nrow(ambiguous_edges) > 0L) {
    p <- p + ggplot2::geom_segment(
      data = ambiguous_edges,
      mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = plot_width),
      inherit.aes = FALSE,
      colour = "grey55",
      linetype = 2,
      alpha = 0.8
    )
  }

  if (nrow(directed_edges) > 0L) {
    p <- p + ggplot2::geom_segment(
      data = directed_edges,
      mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = plot_width),
      inherit.aes = FALSE,
      colour = "black",
      alpha = 0.85,
      arrow = grid::arrow(length = grid::unit(arrow_size, "inches"), type = "closed")
    )
  }

  if (identical(node_colour_by, "cluster")) {
    p <- p + ggplot2::geom_point(
      data = node_df,
      mapping = ggplot2::aes(x = x, y = y, size = plot_size, fill = cluster),
      inherit.aes = FALSE,
      shape = 21,
      colour = "black",
      stroke = 0.5
    ) +
      ggplot2::scale_fill_manual(values = cluster_palette, name = "cluster")
  } else {
    p <- p + ggplot2::geom_point(
      data = node_df,
      mapping = ggplot2::aes(x = x, y = y, size = plot_size, fill = median_time),
      inherit.aes = FALSE,
      shape = 21,
      colour = "black",
      stroke = 0.5
    ) +
      ggplot2::scale_fill_gradientn(
        colours = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"),
        name = "median_time"
      )
  }

  p <- p + ggplot2::geom_text(
    data = node_df,
    mapping = ggplot2::aes(x = x, y = y, label = cluster),
    inherit.aes = FALSE,
    size = vertex_label_size,
    nudge_y = 0.12
  )

  p <- p +
    ggplot2::scale_size_identity() +
    ggplot2::scale_linewidth_identity() +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = main,
      x = if (identical(layout_mode, "embedding")) paste0(dimred, "_1") else "graph_1",
      y = if (identical(layout_mode, "embedding")) paste0(dimred, "_2") else "graph_2"
    )

  p
}


#' Plot TATA trajectories over an embedding
#'
#' Overlay inferred TATA trajectory paths on top of a reduced-dimensional
#' embedding such as UMAP. The plotted paths follow root-to-terminal routes in
#' the cluster graph and can be shown as smoothed curves for a more classical
#' trajectory-style visualization.
#'
#' @param tata_result Result list returned by `run_tata()`.
#' @param dimred Reduced dimension to plot.
#' @param colour_by Optional `colData` column used to color cells. Set to `NULL`
#'   for a grey background.
#' @param mode Display mode for the inferred trajectories. Choices are:
#' - `"full"`: draw one root-to-terminal path for each terminal cluster.
#' - `"combined"`: compress paths with a shared trunk and draw one
#'   representative path per major branch.
#' @param max_paths Optional maximum number of paths to plot after path
#'   selection. `NULL` keeps all available paths.
#' @param point_size Background point size.
#' @param point_alpha Background point alpha.
#' @param smooth_paths Whether to smooth the inferred trajectory paths.
#' @param n_curve_points Number of interpolation points per curve.
#' @param show_centroids Whether to draw cluster centroids.
#' @param show_labels Whether to label cluster centroids.
#' @param curve_colour Colour of the inferred trajectory curves.
#' @param label_size Cluster label size.
#'
#' @return A `ggplot2` object.
#' @export
plot_tata_trajectory_embedding <- function(
    tata_result,
    dimred = "UMAP",
    colour_by = "timepoint",
    mode = c("full", "combined"),
    max_paths = NULL,
    point_size = 0.6,
    point_alpha = 0.75,
    smooth_paths = TRUE,
    n_curve_points = 100,
    show_centroids = TRUE,
    show_labels = TRUE,
    curve_colour = "black",
    label_size = 3) {
  mode <- match.arg(mode)

  sce <- tata_result$sce
  if (!methods::is(sce, "SingleCellExperiment")) {
    stop("`tata_result$sce` must be a SingleCellExperiment.", call. = FALSE)
  }

  if (!is.null(max_paths)) {
    max_paths <- as.integer(max_paths)
    if (length(max_paths) != 1L || is.na(max_paths) || max_paths < 1L) {
      stop("`max_paths` must be NULL or a single positive integer.", call. = FALSE)
    }
  }

  if (!dimred %in% SingleCellExperiment::reducedDimNames(sce)) {
    stop("Reduced dimension `", dimred, "` is not present in `tata_result$sce`.", call. = FALSE)
  }

  cluster_col <- tata_result$parameters$cluster_col
  embedding <- as.matrix(SingleCellExperiment::reducedDim(sce, dimred))
  if (ncol(embedding) < 2L) {
    stop("The requested reduced dimension must contain at least two columns.", call. = FALSE)
  }

  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  meta$Dim1 <- embedding[, 1]
  meta$Dim2 <- embedding[, 2]
  meta$cluster_plot <- meta[[cluster_col]]

  cluster_centroids <- stats::aggregate(
    meta[, c("Dim1", "Dim2")],
    by = list(cluster = as.character(meta$cluster_plot)),
    FUN = stats::median
  )

  cluster_pt <- tata_result$cluster_pseudotime
  if (!is.null(cluster_pt) && nrow(cluster_pt) > 0L) {
    cluster_centroids <- merge(cluster_centroids, cluster_pt, by = "cluster", all.x = TRUE, sort = FALSE)
  } else {
    cluster_centroids$pseudotime <- NA_real_
    cluster_centroids$scaled_pseudotime <- NA_real_
    cluster_centroids$reachable <- TRUE
    cluster_centroids$is_root <- FALSE
  }

  cluster_graph <- tata_result$cluster_graph
  reachable_clusters <- cluster_centroids$cluster[(cluster_centroids$reachable %in% TRUE) | is.na(cluster_centroids$reachable)]
  reachable_clusters <- intersect(reachable_clusters, igraph::V(cluster_graph)$name)

  if (!any(cluster_centroids$is_root, na.rm = TRUE)) {
    root_cluster <- cluster_centroids$cluster[which.min(cluster_centroids$pseudotime)]
  } else {
    root_cluster <- cluster_centroids$cluster[which(cluster_centroids$is_root)[1]]
  }

  subgraph <- igraph::induced_subgraph(cluster_graph, vids = reachable_clusters)
  undirected_subgraph <- igraph::as_undirected(subgraph, mode = "collapse")

  edge_table <- tata_result$cluster_edge_table
  directed_edges <- edge_table[
    edge_table$kept & edge_table$direction != "ambiguous",
    c("from", "to"),
    drop = FALSE
  ]

  cluster_pt_lookup <- stats::setNames(cluster_centroids$pseudotime, cluster_centroids$cluster)

  reachable_order <- cluster_centroids$cluster[
    order(cluster_centroids$pseudotime, decreasing = TRUE, na.last = NA)
  ]
  reachable_order <- intersect(reachable_order, reachable_clusters)

  terminal_clusters <- reachable_order[vapply(
    reachable_order,
    FUN.VALUE = logical(1),
    FUN = function(cluster_name) {
      if (identical(cluster_name, root_cluster)) {
        return(FALSE)
      }

      outgoing <- directed_edges$to[directed_edges$from == cluster_name]
      outgoing <- intersect(outgoing, reachable_clusters)
      if (length(outgoing) == 0L) {
        return(TRUE)
      }

      later_outgoing <- outgoing[cluster_pt_lookup[outgoing] > cluster_pt_lookup[cluster_name]]
      length(later_outgoing) == 0L
    }
  )]

  if (length(terminal_clusters) == 0L) {
    terminal_clusters <- utils::head(setdiff(reachable_order, root_cluster), 3)
  }

  path_info_list <- vector("list", length(terminal_clusters))
  for (i in seq_along(terminal_clusters)) {
    terminal_cluster <- terminal_clusters[i]
    path_vertices <- suppressWarnings(
      igraph::shortest_paths(
        subgraph,
        from = root_cluster,
        to = terminal_cluster,
        mode = "out",
        weights = 1 / pmax(igraph::E(subgraph)$weight, 1e-8)
      )$vpath[[1]]
    )

    if (length(path_vertices) == 0L) {
      path_vertices <- suppressWarnings(
        igraph::shortest_paths(
          undirected_subgraph,
          from = root_cluster,
          to = terminal_cluster,
          mode = "all",
          weights = 1 / pmax(igraph::E(undirected_subgraph)$weight, 1e-8)
        )$vpath[[1]]
      )
    }

    path_clusters <- igraph::V(subgraph)$name[as.integer(path_vertices)]
    path_df <- cluster_centroids[match(path_clusters, cluster_centroids$cluster), c("cluster", "Dim1", "Dim2", "scaled_pseudotime"), drop = FALSE]
    path_df <- path_df[stats::complete.cases(path_df[, c("Dim1", "Dim2")]), , drop = FALSE]
    if (nrow(path_df) < 2L) {
      next
    }

    terminal_pt <- cluster_pt_lookup[terminal_cluster]
    if (is.na(terminal_pt)) {
      terminal_pt <- nrow(path_df)
    }

    path_info_list[[i]] <- list(
      terminal_cluster = terminal_cluster,
      terminal_pseudotime = terminal_pt,
      path_clusters = path_clusters,
      path_df = path_df
    )
  }

  path_info_list <- Filter(Negate(is.null), path_info_list)

  if (identical(mode, "combined") && length(path_info_list) > 1L) {
    min_path_len <- min(vapply(path_info_list, function(x) length(x$path_clusters), integer(1)))
    common_prefix_len <- 0L
    for (pos in seq_len(min_path_len)) {
      prefix_values <- unique(vapply(path_info_list, function(x) x$path_clusters[pos], character(1)))
      if (length(prefix_values) == 1L) {
        common_prefix_len <- pos
      } else {
        break
      }
    }

    branch_key <- vapply(
      path_info_list,
      FUN.VALUE = character(1),
      FUN = function(path_info) {
        if (length(path_info$path_clusters) > common_prefix_len) {
          path_info$path_clusters[common_prefix_len + 1L]
        } else {
          utils::tail(path_info$path_clusters, 1)
        }
      }
    )

    grouped_paths <- split(seq_along(path_info_list), branch_key)
    keep_idx <- vapply(
      grouped_paths,
      FUN.VALUE = integer(1),
      FUN = function(idx) {
        candidate_pt <- vapply(path_info_list[idx], function(x) x$terminal_pseudotime, numeric(1))
        idx[which.max(candidate_pt)]
      }
    )
    path_info_list <- path_info_list[unname(keep_idx)]
  }

  if (!is.null(max_paths) && length(path_info_list) > max_paths) {
    path_rank <- order(
      vapply(path_info_list, function(x) x$terminal_pseudotime, numeric(1)),
      decreasing = TRUE
    )
    path_info_list <- path_info_list[path_rank[seq_len(max_paths)]]
  }

  curve_list <- vector("list", length(path_info_list))
  for (i in seq_along(path_info_list)) {
    path_df <- path_info_list[[i]]$path_df

    if (isTRUE(smooth_paths) && nrow(path_df) >= 4L) {
      spline_index <- seq_len(nrow(path_df))
      x_spline <- stats::smooth.spline(spline_index, path_df$Dim1, spar = 0.5)
      y_spline <- stats::smooth.spline(spline_index, path_df$Dim2, spar = 0.5)
      new_index <- seq(min(spline_index), max(spline_index), length.out = n_curve_points)
      curve_df <- data.frame(
        path_id = paste0("path_", i),
        terminal_cluster = path_info_list[[i]]$terminal_cluster,
        Dim1 = stats::predict(x_spline, x = new_index)$y,
        Dim2 = stats::predict(y_spline, x = new_index)$y,
        stringsAsFactors = FALSE
      )
    } else {
      curve_df <- data.frame(
        path_id = paste0("path_", i),
        terminal_cluster = path_info_list[[i]]$terminal_cluster,
        Dim1 = path_df$Dim1,
        Dim2 = path_df$Dim2,
        stringsAsFactors = FALSE
      )
    }

    curve_list[[i]] <- curve_df
  }

  curve_df <- do.call(rbind, Filter(Negate(is.null), curve_list))

  if (!is.null(colour_by) && colour_by %in% colnames(meta)) {
    meta$colour_value <- meta[[colour_by]]
    discrete_colour <- !is.numeric(meta$colour_value)
  } else {
    meta$colour_value <- "cells"
    discrete_colour <- TRUE
  }

  p <- ggplot2::ggplot(meta, ggplot2::aes(x = Dim1, y = Dim2))
  if (!is.null(colour_by) && colour_by %in% colnames(meta)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(colour = colour_value),
      size = point_size,
      alpha = point_alpha
    )
  } else {
    p <- p + ggplot2::geom_point(
      colour = "grey80",
      size = point_size,
      alpha = point_alpha
    )
  }

  if (!is.null(curve_df) && nrow(curve_df) > 0L) {
    p <- p + ggplot2::geom_path(
      data = curve_df,
      mapping = ggplot2::aes(x = Dim1, y = Dim2, group = path_id),
      inherit.aes = FALSE,
      colour = curve_colour,
      linewidth = 1.1,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.15, "inches"), type = "closed")
    )
  }

  if (isTRUE(show_centroids)) {
    p <- p + ggplot2::geom_point(
      data = cluster_centroids,
      mapping = ggplot2::aes(x = Dim1, y = Dim2),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.6,
      fill = "white",
      colour = "black",
      stroke = 0.6
    )
  }

  if (isTRUE(show_labels)) {
    p <- p + ggplot2::geom_text(
      data = cluster_centroids,
      mapping = ggplot2::aes(x = Dim1, y = Dim2, label = cluster),
      inherit.aes = FALSE,
      size = label_size,
      nudge_y = 0.15
    )
  }

  p <- p +
    ggplot2::theme_classic() +
    ggplot2::labs(x = paste0(dimred, "_1"), y = paste0(dimred, "_2"), colour = colour_by)

  if (!is.null(colour_by) && colour_by %in% colnames(meta) && !discrete_colour) {
    p <- p + ggplot2::scale_colour_viridis_c(option = "plasma")
  }

  p
}
