# Column-wise tied ranks of a matrix

Ranks the entries of each column independently using the average rank
for ties, and returns the per-column tie group sizes needed for the
Wilcoxon variance correction. Used internally by
[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
(on a transposed input, so that rows become observations) but exposed as
a fast standalone ranking primitive for sparse and dense numeric
matrices.

## Usage

``` r
rank_matrix(X)

# S3 method for class 'dgCMatrix'
rank_matrix(X)

# S3 method for class 'matrix'
rank_matrix(X)
```

## Arguments

- X:

  Numeric matrix or `dgCMatrix`.

## Value

List with two elements:

- `X_ranked` - matrix with the same shape as `X` containing per-column
  tied ranks.

- `ties` - list of integer vectors, one per column, giving the sizes of
  all tie groups encountered in that column. Used by the Wilcoxon
  statistic to correct for ties.

## See also

[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)

## Examples

``` r
set.seed(42)
exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
                dimnames = list(paste0("G", 1:25), NULL))
rank_res <- rank_matrix(exprs)
```
