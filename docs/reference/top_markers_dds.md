# Top markers per group from pseudobulk DESeq2 results

Filters and ranks the long-form output of
[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md)
to give the most distinguishing features per group. The filtering
arguments combine multiplicatively, then the top `n` features per group
are kept by descending Wald statistic and pivoted into wide form.
Counterpart to
[`top_markers()`](https://immunogenomics.github.io/presto/reference/top_markers.md)
for Wilcoxon-based results.

## Usage

``` r
top_markers_dds(res, n = 10, pval_max = 1, padj_max = 1, lfc_min = 1)
```

## Arguments

- res:

  Long-form DESeq2 results from
  [`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md).
  Must have columns `group`, `feature`, `pvalue`, `padj`,
  `log2FoldChange`, and `stat`.

- n:

  Number of top features to return per group. Default `10`.

- pval_max:

  Filter features with raw `pvalue > pval_max`. Default `1` (no filter).

- padj_max:

  Filter features with adjusted `padj > padj_max`. Default `1` (no
  filter).

- lfc_min:

  Filter features with `log2FoldChange < lfc_min`. Default `1` keeps
  only upregulated features; set to `0` to include downregulated as
  well, or `-Inf` to disable.

## Value

tibble in wide form: a `rank` column (1..`n`) and one column per group
containing the gene identifier of the top-ranked feature at that rank.
Cells are `NA` for groups that have fewer than `n` features passing the
filters.

## See also

[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md),
[`top_markers()`](https://immunogenomics.github.io/presto/reference/top_markers.md)
