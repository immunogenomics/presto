# Get top n markers from pseudobulk DESeq2

Useful summary of the most distinguishing features in each group.

## Usage

``` r
top_markers_dds(res, n = 10, pval_max = 1, padj_max = 1, lfc_min = 1)
```

## Arguments

- res:

  table returned by pseudobulk_deseq2() function.

- n:

  number of markers to find for each.

- pval_max:

  filter features with pval \> pval_max.

- padj_max:

  filter features with padj \> padj_max.

- lfc_min:

  filter features with log2FoldChange \< lfc_min

  @return table with the top n markers for each cluster.
