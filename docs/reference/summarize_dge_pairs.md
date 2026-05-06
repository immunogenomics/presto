# Summarize directional pairwise DGE results

Collapses the long-form output of `pseudobulk_deseq2(mode = "pairwise")`
to a single row per `(group, gene)` by selecting one comparison per
gene-group pair. With `mode = "min"` it keeps the worst comparison (most
conservative — the gene must beat *every* other group to have a high
statistic). With `mode = "max"` it keeps the best (most permissive — the
gene need only beat one).

## Usage

``` r
summarize_dge_pairs(dge_res, mode = c("min", "max")[1])
```

## Arguments

- dge_res:

  data.frame of pairwise results from
  [`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md)
  with `mode = "pairwise"`. Must have columns `group1`, `group2`,
  `feature`, and `stat`.

- mode:

  `"min"` (default) keeps the comparison with the smallest `stat` for
  each (group, gene); `"max"` keeps the largest.

## Value

data.frame with one row per (group, gene), sorted by `stat` descending
within each group. The `group2` column is dropped and `group1` is
renamed to `group`.

## See also

[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md),
[`pseudobulk_pairwise()`](https://immunogenomics.github.io/presto/reference/pseudobulk_pairwise.md)
