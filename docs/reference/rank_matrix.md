# rank_matrix

Utility function to rank columns of matrix

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

  feature by observation matrix.

## Value

List with 2 items

- X_ranked - matrix of entry ranks

- ties - list of tied group sizes

## Examples

``` r

data(exprs)
rank_res <- rank_matrix(exprs)
```
