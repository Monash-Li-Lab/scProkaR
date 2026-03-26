# Plot TATA trajectories over an embedding

Overlay inferred TATA trajectory paths on top of a reduced-dimensional
embedding such as UMAP. The plotted paths follow root-to-terminal routes
in the cluster graph and can be shown as smoothed curves for a more
classical trajectory-style visualization.

## Usage

``` r
plot_tata_trajectory_embedding(
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
  label_size = 3
)
```

## Arguments

- tata_result:

  Result list returned by
  [`run_tata()`](https://monash-li-lab.github.io/scProka/reference/run_tata.md).

- dimred:

  Reduced dimension to plot.

- colour_by:

  Optional `colData` column used to color cells. Set to `NULL` for a
  grey background.

- mode:

  Display mode for the inferred trajectories. Choices are:

  - `"full"`: draw one root-to-terminal path for each terminal cluster.

  - `"combined"`: compress paths with a shared trunk and draw one
    representative path per major branch.

- max_paths:

  Optional maximum number of paths to plot after path selection. `NULL`
  keeps all available paths.

- point_size:

  Background point size.

- point_alpha:

  Background point alpha.

- smooth_paths:

  Whether to smooth the inferred trajectory paths.

- n_curve_points:

  Number of interpolation points per curve.

- show_centroids:

  Whether to draw cluster centroids.

- show_labels:

  Whether to label cluster centroids.

- curve_colour:

  Colour of the inferred trajectory curves.

- label_size:

  Cluster label size.

## Value

A `ggplot2` object.
