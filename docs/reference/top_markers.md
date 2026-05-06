# Get top n markers from wilcoxauc

Useful summary of the most distinguishing features in each group.

## Usage

``` r
top_markers(
  res,
  n = 10,
  auc_min = 0,
  pval_max = 1,
  padj_max = 1,
  pct_in_min = 0,
  pct_out_max = 100
)
```

## Arguments

- res:

  table returned by wilcoxauc() function.

- n:

  number of markers to find for each.

- auc_min:

  filter features with auc \< auc_min.

- pval_max:

  filter features with pval \> pval_max.

- padj_max:

  filter features with padj \> padj_max.

- pct_in_min:

  Minimum percent (0-100) of observations with non-zero entries in
  group.

- pct_out_max:

  Maximum percent (0-100) of observations with non-zero entries out of
  group.

## Value

table with the top n markers for each cluster.

## Examples

``` r

data(exprs)
data(y)

## first, run wilcoxauc
res <- wilcoxauc(exprs, y)

## top 10 markers for each group
## filter for nominally significant (p<0.05) and over-expressed (auc>0.5)
top_markers(res, 10, auc_min = 0.5, pval_max = 0.05)
#> # A tibble: 9 × 4
#>    rank A     B     C    
#>   <int> <chr> <chr> <chr>
#> 1     1 G4    G20   G1   
#> 2     2 NA    G5    G21  
#> 3     3 NA    G15   G6   
#> 4     4 NA    G19   G16  
#> 5     5 NA    G25   G11  
#> 6     6 NA    G10   G7   
#> 7     7 NA    NA    G2   
#> 8     8 NA    NA    G17  
#> 9     9 NA    NA    G12  
```
