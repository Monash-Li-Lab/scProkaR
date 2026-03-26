# Orient pseudotime to match increasing real time

Flip an inferred pseudotime vector when needed so that it best agrees
with increasing experimental time. This is useful when benchmarking
methods whose pseudotime direction is arbitrary up to reversal.

## Usage

``` r
orient_pseudotime_to_time(
  pseudotime,
  timepoint,
  method = c("spearman", "kendall", "pearson")
)
```

## Arguments

- pseudotime:

  Numeric vector of inferred pseudotime values.

- timepoint:

  Experimental time labels for the same cells.

- method:

  Correlation used to choose the orientation. Choices are:

  - `"spearman"`: rank-based and the default choice for monotone
    trajectories.

  - `"kendall"`: more conservative rank correlation.

  - `"pearson"`: linear correlation on the scaled values.

## Value

A numeric vector scaled to `[0, 1]` and oriented so that larger values
agree as well as possible with later experimental time.
