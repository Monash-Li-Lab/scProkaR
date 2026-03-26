# Merge multiple SCProkaR or SingleCellExperiment objects

`MergeBacObjects()` combines multiple per-sample `SingleCellExperiment`
objects without relying on external
[`cbind()`](https://rdrr.io/r/base/cbind.html) method dispatch. This is
useful when merging sample-level objects created by
[`CreateBacObject()`](https://monash-li-lab.github.io/scProka/reference/CreateBacObject.md)
before batch integration.

## Usage

``` r
MergeBacObjects(..., objects = NULL, gene_mode = c("intersect", "union"))
```

## Arguments

- ...:

  Individual `SingleCellExperiment` objects to merge.

- objects:

  Optional list of `SingleCellExperiment` objects.

- gene_mode:

  Either `"intersect"` or `"union"`.

## Value

A merged `SingleCellExperiment`.

## Details

By default the function keeps the intersection of genes across inputs,
which is the safest choice before integration.
