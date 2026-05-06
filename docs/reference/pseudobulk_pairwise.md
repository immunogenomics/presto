# Pseudobulk DESeq2: pairwise contrasts

For each ordered pair of contrast levels in `vals_test`, fits a DESeq2
model on just those two groups of pseudobulks and returns the
foreground-vs-background coefficient. This is more conservative than
one-vs-all because a marker has to differentiate the foreground from
*every* other group, not just the average background. Cost grows as
O(N^2) in the number of levels, so subset to the levels of interest
first.

## Usage

``` r
pseudobulk_pairwise(
  dge_formula,
  counts_df,
  meta_data,
  contrast_var,
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

- contrast_var:

  Name of the contrast column in `meta_data` (the first term of
  `dge_formula`).

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

data.frame of DESeq2 results with columns `group1`, `group2`, `feature`,
`baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`. Each
row reports the test where pseudobulks of `group1` are the foreground
and pseudobulks of `group2` are the background.

## Details

Most users should call
[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md)
with `mode = "pairwise"` rather than this function directly. The
companion
[`summarize_dge_pairs()`](https://immunogenomics.github.io/presto/reference/summarize_dge_pairs.md)
collapses the directional pairs to a single row per (gene, group).

## See also

[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md),
[`pseudobulk_one_vs_all()`](https://immunogenomics.github.io/presto/reference/pseudobulk_one_vs_all.md),
[`pseudobulk_within()`](https://immunogenomics.github.io/presto/reference/pseudobulk_within.md),
[`summarize_dge_pairs()`](https://immunogenomics.github.io/presto/reference/summarize_dge_pairs.md)
