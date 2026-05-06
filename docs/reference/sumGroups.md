# Group-wise sum of a matrix along one axis

For each unique value of the grouping vector `y`, sums the corresponding
rows (or columns) of `X`. Used internally by
[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
and
[`collapse_counts()`](https://immunogenomics.github.io/presto/reference/collapse_counts.md),
but exposed as a fast group-wise reduction primitive that works on both
dense matrices and `dgCMatrix` sparse inputs.

## Usage

``` r
sumGroups(X, y, MARGIN = 2)

# S3 method for class 'dgCMatrix'
sumGroups(X, y, MARGIN = 2)

# S3 method for class 'matrix'
sumGroups(X, y, MARGIN = 2)
```

## Arguments

- X:

  Numeric matrix or `dgCMatrix`.

- y:

  Group label vector. Coerced to integer factor codes.

- MARGIN:

  Whether observations are along rows or columns of `X`. `MARGIN = 2`
  (default): observations are rows (`length(y) == nrow(X)`); rows are
  summed within each group. `MARGIN = 1`: observations are columns
  (`length(y) == ncol(X)`); columns are summed within each group.

## Value

Numeric matrix of shape `n_groups x n_features`. Row order matches the
integer order of `factor(y)`.

## See also

[`nnzeroGroups()`](https://immunogenomics.github.io/presto/reference/nnzeroGroups.md),
[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)

## Examples

``` r
set.seed(42)
exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
                dimnames = list(paste0("G", 1:25), NULL))
y <- rep(c("A", "B", "C"), each = 50)
sumGroups_res <- sumGroups(exprs, y, 1)
#> Warning: NAs introduced by coercion
sumGroups_res <- sumGroups(t(exprs), y, 2)
#> Warning: NAs introduced by coercion
```
