# Run UMAP on an existing embedding

Computes a two-dimensional UMAP from a stored reduced dimension, which
is useful for inspecting integrated embeddings returned by
[`IntegrateBacData()`](https://monash-li-lab.github.io/scProka/reference/IntegrateBacData.md)
or
[`RegisterIntegrationEmbedding()`](https://monash-li-lab.github.io/scProka/reference/RegisterIntegrationEmbedding.md).

## Usage

``` r
RunIntegratedUMAP(
  sce,
  reduction,
  umap_name = NULL,
  n_neighbors = 30,
  min_dist = 0.3,
  metric = "cosine",
  seed = 1,
  ...
)
```

## Arguments

- sce:

  A `SingleCellExperiment`.

- reduction:

  Existing reduced-dimension name to use as UMAP input.

- umap_name:

  Optional name for the output UMAP reduced dimension.

- n_neighbors:

  Number of neighbours passed to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).

- min_dist:

  UMAP `min_dist`.

- metric:

  Distance metric passed to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).

- seed:

  Random seed for reproducibility.

- ...:

  Additional arguments passed to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).

## Value

A `SingleCellExperiment` with the computed UMAP stored in
`reducedDims(sce)`.
