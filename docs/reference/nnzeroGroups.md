# Group-wise non-zero counts of a matrix along one axis

For each unique value of the grouping vector `y`, counts the number of
non-zero entries among the corresponding rows (or columns) of `X`. Used
internally by
[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
to compute the percent-expressed columns (`pct_in`, `pct_out`), but
exposed as a fast group-wise reduction primitive for both dense and
`dgCMatrix` inputs.

## Usage

``` r
nnzeroGroups(X, y, MARGIN = 2)

# S3 method for class 'dgCMatrix'
nnzeroGroups(X, y, MARGIN = 2)

# S3 method for class 'matrix'
nnzeroGroups(X, y, MARGIN = 2)
```

## Arguments

- X:

  Numeric matrix or `dgCMatrix`.

- y:

  Group label vector. Coerced to integer factor codes.

- MARGIN:

  Whether observations are along rows or columns of `X`. `MARGIN = 2`
  (default): observations are rows (`length(y) == nrow(X)`).
  `MARGIN = 1`: observations are columns (`length(y) == ncol(X)`).

## Value

Integer matrix of shape `n_groups x n_features`, where entry `(g, j)` is
the number of observations in group `g` for which feature `j` is
non-zero.

## See also

[`sumGroups()`](https://immunogenomics.github.io/presto/reference/sumGroups.md),
[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)

## Examples

``` r
set.seed(42)
exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
                dimnames = list(paste0("G", 1:25), NULL))
y <- rep(c("A", "B", "C"), each = 50)
nnz_res <- nnzeroGroups(exprs, y, 1)
#> Warning: NAs introduced by coercion
nnz_res <- nnzeroGroups(t(exprs), y, 2)
#> Warning: NAs introduced by coercion
```
