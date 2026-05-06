# sumGroups

Utility function to sum over group labels

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
sumGroups_res <- sumGroups(exprs, y, 1)
sumGroups_res <- sumGroups(t(exprs), y, 2)
```
