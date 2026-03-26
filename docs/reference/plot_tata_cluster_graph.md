# Plot the cluster-level TATA graph

Plot the directed cluster-level TATA graph. The graph can be drawn in
its own abstract layout or positioned on a chosen reduced dimension so
that the graph coordinates match the embedding shown to the user.

## Usage

``` r
plot_tata_cluster_graph(
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
  arrow_size = 0.15
)
```

## Arguments

- tata_result:

  Result list returned by
  [`run_tata()`](https://monash-li-lab.github.io/scProka/reference/run_tata.md).

- dimred:

  Optional reduced dimension used to position cluster nodes. When
  supplied together with `layout_mode = "embedding"`, cluster centroids
  are computed from this embedding.

- layout_mode:

  Layout mode for the graph. Choices are:

  - `"graph"`: use the stored abstract graph coordinates, or a
    force-directed layout if they are unavailable.

  - `"embedding"`: place cluster nodes at their centroids in the
    selected reduced dimension.

- cluster_col:

  Optional cluster column used when `layout_mode = "embedding"`. By
  default this is taken from `tata_result$parameters`.

- show_cells:

  Logical indicating whether to draw cells in the background when
  `layout_mode = "embedding"`.

- colour_by:

  Optional `colData` column used to color background cells. If `NULL`
  and `show_cells = TRUE`, the cluster column is used when available.

- node_colour_by:

  Node coloring scheme. Choices are:

  - `"median_time"`: color nodes by cluster median time.

  - `"cluster"`: color nodes by cluster identity, which is useful for
    embedding-aligned graph plots.

- main:

  Plot title.

- vertex_label_size:

  Label size.

- edge_width_scale:

  Relative edge width scale.

- point_size:

  Background point size when cells are drawn.

- point_alpha:

  Background point alpha when cells are drawn.

- arrow_size:

  Arrow size for directed edges.

## Value

A `ggplot2` object.
