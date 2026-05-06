# Pseudobulk DESeq2: within-group contrasts

Splits the pseudobulks by `split_var` and, within each split level, fits
a DESeq2 model whose contrast variable is the second term of
`dge_formula`. Used to test for an effect (e.g. case vs. control)
restricted to one cluster at a time, where pooling across clusters would
mix biology and batch.

## Usage

``` r
pseudobulk_within(
  dge_formula,
  counts_df,
  meta_data,
  split_var,
  vals_test,
  verbose,
  min_counts_per_sample,
  present_in_min_samples
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

- split_var:

  Name of the column in `meta_data` to split on (the first term of
  `dge_formula`). Each unique level of `split_var` produces an
  independent DESeq2 fit.

- vals_test:

  Character vector of contrast levels to test. If `NULL` (default),
  every level of the contrast variable is tested.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

- min_counts_per_sample:

  Minimum count per pseudobulk for a gene to be considered expressed in
  that pseudobulk. Default `10`.

- present_in_min_samples:

  Minimum number of pseudobulks in which a gene must reach
  `min_counts_per_sample` to be retained. Default `5`.

## Value

data.frame of DESeq2 results with columns `group` (the value of
`split_var`), `feature`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`,
`pvalue`, `padj`.

## Details

Two-level contrasts (factor or character) yield the standard
`<var>_<level2>_vs_<level1>` Wald coefficient. Three or more levels are
integer-encoded and the fit returns an ordinal trend. Most users should
call
[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md)
with `mode = "within"` rather than this function directly. See the
`pseudobulk` vignette for a worked example using case-vs-control DGE
within each cell cluster.

## See also

[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md),
[`pseudobulk_one_vs_all()`](https://immunogenomics.github.io/presto/reference/pseudobulk_one_vs_all.md),
[`pseudobulk_pairwise()`](https://immunogenomics.github.io/presto/reference/pseudobulk_pairwise.md)
