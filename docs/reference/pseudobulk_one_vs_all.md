# Pseudobulk DESeq2: one-vs-all contrasts

For each level of `contrast_var`, fits a DESeq2 model where that level
is the foreground and every other level is pooled as the background.
Useful for marker-style differential expression: the return is one set
of `(log2FoldChange, padj)` per gene per group, giving a quick view of
what's elevated in each group relative to the rest.

## Usage

``` r
pseudobulk_one_vs_all(
  dge_formula,
  counts_df,
  meta_data,
  contrast_var,
  vals_test,
  collapse_background,
  verbose
)
```

## Arguments

- dge_formula:

  One-sided formula such as `~cluster + donor`. The first term is
  treated as the contrast variable in `"one_vs_all"` and `"pairwise"`
  modes, or as the split variable in `"within"` mode (in which case the
  second term becomes the contrast variable). Additional terms are kept
  as covariates in the DESeq2 design.

- counts_df:

  Feature-by-pseudobulk integer count matrix. Rows are features; columns
  must align with rows of `meta_data`.

- meta_data:

  data.frame of pseudobulk metadata. One row per pseudobulk; should
  contain only the variables used in `dge_formula`.

- contrast_var:

  Name of the contrast column in `meta_data` (the first term of
  `dge_formula`).

- vals_test:

  Character vector of contrast levels to test. If `NULL` (default),
  every level of the contrast variable is tested.

- collapse_background:

  Used only when `mode = "one_vs_all"`. If `TRUE`, background
  pseudobulks are collapsed across the contrast variable before each
  DESeq2 fit, which helps when donor representation is unbalanced across
  clusters. Default `TRUE`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

data.frame of DESeq2 results with columns `group`, `feature`,
`baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`. Sorted
by `stat` descending within each group.

## Details

Most users should call
[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md)
with `mode = "one_vs_all"` rather than this function directly; the
wrapper handles formula parsing, gene-count filtering, and dispatch. See
the `pseudobulk` vignette for a worked example.

## See also

[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md),
[`pseudobulk_pairwise()`](https://immunogenomics.github.io/presto/reference/pseudobulk_pairwise.md),
[`pseudobulk_within()`](https://immunogenomics.github.io/presto/reference/pseudobulk_within.md),
[`top_markers_dds()`](https://immunogenomics.github.io/presto/reference/top_markers_dds.md)
