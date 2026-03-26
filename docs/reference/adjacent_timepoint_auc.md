# Compute adjacent-timepoint AUC for a pseudotime ordering

For each adjacent pair of experimental timepoints, compute the
probability that a cell from the later timepoint has larger pseudotime
than a cell from the earlier timepoint. This is equivalent to the
Mann-Whitney AUC.

## Usage

``` r
adjacent_timepoint_auc(
  pseudotime,
  timepoint,
  group = NULL,
  average = c("weighted", "macro"),
  min_cells = 10L
)
```

## Arguments

- pseudotime:

  Numeric vector of inferred pseudotime values.

- timepoint:

  Experimental time labels for the same cells.

- group:

  Optional grouping variable, for example treatment or condition. If
  provided, adjacent-timepoint AUC is computed within each group before
  averaging.

- average:

  How to summarize the per-pair AUC values. Choices are:

  - `"weighted"`: weight each comparison by the number of earlier-later
    cell pairs.

  - `"macro"`: give every valid comparison equal weight.

- min_cells:

  Minimum number of cells required in each timepoint group for a
  comparison to be included.

## Value

A list with:

- `summary`: a one-row `data.frame` containing the overall
  adjacent-timepoint AUC and the number of valid comparisons.

- `pair_table`: a `data.frame` with AUC values for each adjacent
  timepoint comparison.
