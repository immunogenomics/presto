# Collapse a single-cell count matrix into pseudobulks

Sums (or averages) the columns of a feature-by-cell count matrix
according to one or more cell-metadata columns. The result is a
feature-by-pseudobulk matrix where each column pools all cells that
share the same combination of metadata values (e.g. all cells from a
given donor in a given cluster). The resulting matrix is suitable for
bulk-RNA-seq tools such as DESeq2, edgeR, or limma. See the `pseudobulk`
vignette for an end-to-end walkthrough.

## Usage

``` r
collapse_counts(
  counts_mat,
  meta_data,
  varnames,
  min_cells_per_group = 0,
  keep_n = FALSE,
  how = c("sum", "mean")[1]
)
```

## Arguments

- counts_mat:

  Counts matrix. Rows are features (genes), columns are cells. Sparse
  (`dgCMatrix`) or dense.

- meta_data:

  data.frame of cell metadata. Must have one row per column of
  `counts_mat`.

- varnames:

  Character vector of column names in `meta_data` that together define
  the pseudobulk grouping. Each unique combination of values across
  these columns becomes one pseudobulk sample.

- min_cells_per_group:

  Drop pseudobulks containing fewer than this many cells. Default `0`
  (keep all).

- keep_n:

  If `TRUE`, retain the per-pseudobulk cell count as column `N` in the
  returned `meta_data`. Default `FALSE`.

- how:

  `"sum"` (default) sums counts across cells in each pseudobulk;
  `"mean"` divides by `N` to give the mean expression per cell. `DESeq2`
  and other count-based regressions expect sums.

## Value

A list with two elements:

- `counts_mat` - feature-by-pseudobulk numeric matrix.

- `meta_data` - data.frame with one row per pseudobulk containing the
  columns named in `varnames` (and `N` if `keep_n = TRUE`).

## See also

[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md),
[`compute_hash()`](https://immunogenomics.github.io/presto/reference/compute_hash.md)

## Examples

``` r
m <- matrix(sample.int(8, 100*500, replace=TRUE), nrow=100, ncol=500)
rownames(m) <- paste0("G", 1:100)
colnames(m) <- paste0("C", 1:500)
md1 <- sample(c("a", "b"), 500, replace=TRUE)
md2 <- sample(c("c", "d"), 500, replace=TRUE)
df <- data.frame(md1, md2)
data_collapsed <- collapse_counts(m, df, c("md1", "md2"))
head(data_collapsed$counts_mat)
#>    sample_1 sample_0 sample_2 sample_3
#> G1      627      523      550      573
#> G2      663      498      559      554
#> G3      574      517      575      502
#> G4      605      511      551      520
#> G5      639      548      577      546
#> G6      582      554      596      521
head(data_collapsed$meta_data)
#>   md1 md2
#> 1   b   c
#> 2   a   c
#> 3   a   d
#> 4   b   d
```
