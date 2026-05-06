# nnzeroGroups

Utility function to compute number of zeros-per-feature within group

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

  matrix

- y:

  group labels

- MARGIN:

  whether observations are rows (=2) or columns (=1)

## Value

Matrix of groups by features

## Examples

``` r

data(exprs)
data(y)
nnz_res <- nnzeroGroups(exprs, y, 1)
nnz_res <- nnzeroGroups(t(exprs), y, 2)
```
