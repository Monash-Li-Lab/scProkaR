# Measure branch-emergence timing from pseudotime

Compare the known order or timing of branch emergence with the inferred
pseudotime of sufficiently committed cells in each branch.

## Usage

``` r
branch_onset_metrics(
  pseudotime,
  branch,
  branch_activation,
  true_onset,
  activation_threshold = 0.25
)
```

## Arguments

- pseudotime:

  Numeric vector of inferred pseudotime values.

- branch:

  Branch or condition label for each cell.

- branch_activation:

  Numeric activation score indicating how committed a cell is to its
  branch.

- true_onset:

  Named numeric vector giving the known onset time of each branch. The
  names must match a subset of the values in `branch`.

- activation_threshold:

  Minimum activation value required for a cell to be treated as
  branch-committed when estimating the branch onset.

## Value

A list with:

- `summary`: a one-row `data.frame` containing Kendall tau, Spearman
  rho, normalized branch-onset mean absolute error, a higher-is-better
  branch onset order score, a higher-is-better branch onset accuracy
  score, and the number of branches used.

- `onset_table`: a `data.frame` with the true and inferred onset values
  for each branch.
