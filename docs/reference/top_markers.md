# Top markers per group from wilcoxauc results

Filters and ranks the long-form output of
[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
to give the most distinguishing features per group. The filter arguments
combine multiplicatively, then the top `n` features per group are kept
by descending `auc` and pivoted into wide form. Counterpart to
[`top_markers_dds()`](https://immunogenomics.github.io/presto/reference/top_markers_dds.md)
for DESeq2-based pseudobulk results.

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

  Long-form results table from
  [`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md).

- n:

  Number of top markers to return per group. Default `10`.

- auc_min:

  Drop features with `auc < auc_min`. Default `0` (no filter); set to
  `0.5` to keep only features that are positive markers (more highly
  expressed in-group than out).

- pval_max:

  Drop features with raw `pval > pval_max`. Default `1`.

- padj_max:

  Drop features with adjusted `padj > padj_max`. Default `1`.

- pct_in_min:

  Minimum percent (0-100) of in-group observations with non-zero feature
  value. Default `0`.

- pct_out_max:

  Maximum percent (0-100) of out-of-group observations with non-zero
  feature value. Default `100`.

## Value

tibble in wide form: a `rank` column (1..`n`) and one column per group
containing the feature name of the top-ranked marker at that rank. Cells
are `NA` for groups with fewer than `n` features that pass the filters.

## See also

[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md),
[`top_markers_dds()`](https://immunogenomics.github.io/presto/reference/top_markers_dds.md)

## Examples

``` r
set.seed(42)
exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
                dimnames = list(paste0("G", 1:25), NULL))
y <- rep(c("A", "B", "C"), each = 50)

res <- wilcoxauc(exprs, y)

## top 10 markers per group, restricted to nominally significant,
## up-regulated features (auc > 0.5 means in-group > out-of-group).
top_markers(res, n = 10, auc_min = 0.5, pval_max = 0.05)
#> # A tibble: 2 × 2
#>    rank C    
#>   <int> <chr>
#> 1     1 G25  
#> 2     2 G5   
```
