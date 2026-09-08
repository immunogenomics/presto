# Pseudobulk differential expression with DESeq2

Runs [`DESeq2::DESeq()`](https://rdrr.io/pkg/DESeq2/man/DESeq.html) on a
feature-by-pseudobulk count matrix to identify genes whose expression
differs between groups of pseudobulks. Three test designs are supported,
selected by `mode`:

- `"one_vs_all"` (default) - test each level of the contrast variable
  against the union of all others. Useful for marker discovery.
  Dispatched to
  [`pseudobulk_one_vs_all()`](https://immunogenomics.github.io/presto/reference/pseudobulk_one_vs_all.md).

- `"pairwise"` - test every ordered pair of levels separately. Useful
  for high-confidence markers; combine with
  [`summarize_dge_pairs()`](https://immunogenomics.github.io/presto/reference/summarize_dge_pairs.md)
  to keep the most conservative pair. Dispatched to
  [`pseudobulk_pairwise()`](https://immunogenomics.github.io/presto/reference/pseudobulk_pairwise.md).

- `"within"` - split on the first term of `dge_formula` and test the
  second term within each split level. Useful for condition /
  case-vs-control comparisons restricted to one cluster at a time.
  Dispatched to
  [`pseudobulk_within()`](https://immunogenomics.github.io/presto/reference/pseudobulk_within.md).

Pseudobulk inputs are typically produced by
[`collapse_counts()`](https://immunogenomics.github.io/presto/reference/collapse_counts.md).
Genes with low counts across pseudobulks are filtered out before fitting
(controlled by `min_counts_per_sample` and `present_in_min_samples`).
See the `pseudobulk` vignette for an end-to-end walkthrough on a real
single-cell dataset.

## Usage

``` r
pseudobulk_deseq2(
  dge_formula,
  meta_data,
  counts_df,
  verbose = TRUE,
  min_counts_per_sample = 10,
  present_in_min_samples = 5,
  collapse_background = TRUE,
  vals_test = NULL,
  mode = c("one_vs_all", "pairwise", "within")[1]
)
```

## Arguments

- dge_formula:

  One-sided formula such as `~cluster + donor`. The first term is
  treated as the contrast variable in `"one_vs_all"` and `"pairwise"`
  modes, or as the split variable in `"within"` mode (in which case the
  second term becomes the contrast variable). Additional terms are kept
  as covariates in the DESeq2 design.

- meta_data:

  data.frame of pseudobulk metadata. One row per pseudobulk; should
  contain only the variables used in `dge_formula`.

- counts_df:

  Feature-by-pseudobulk integer count matrix. Rows are features; columns
  must align with rows of `meta_data`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

- min_counts_per_sample:

  Minimum count per pseudobulk for a gene to be considered expressed in
  that pseudobulk. Default `10`.

- present_in_min_samples:

  Minimum number of pseudobulks in which a gene must reach
  `min_counts_per_sample` to be retained. Default `5`.

- collapse_background:

  Used only when `mode = "one_vs_all"`. If `TRUE`, background
  pseudobulks are collapsed across the contrast variable before each
  DESeq2 fit, which helps when donor representation is unbalanced across
  clusters. Default `TRUE`.

- vals_test:

  Character vector of contrast levels to test. If `NULL` (default),
  every level of the contrast variable is tested.

- mode:

  One of `"one_vs_all"` (default), `"pairwise"`, or `"within"`. See
  description.

## Value

A long-form data.frame of DESeq2 results. Columns include the group
identifier(s) (`group`, or `group1` / `group2` in `pairwise` mode),
`feature`, and the DESeq2 columns `baseMean`, `log2FoldChange`, `lfcSE`,
`stat`, `pvalue`, `padj`.

## See also

[`collapse_counts()`](https://immunogenomics.github.io/presto/reference/collapse_counts.md),
[`top_markers_dds()`](https://immunogenomics.github.io/presto/reference/top_markers_dds.md),
[`summarize_dge_pairs()`](https://immunogenomics.github.io/presto/reference/summarize_dge_pairs.md)

## Examples

``` r
# \donttest{
if (requireNamespace("DESeq2", quietly = TRUE)) {
    ## 100 genes x 500 cells from 2 clusters across 6 donors
    m <- matrix(sample.int(8, 100 * 500, replace = TRUE), nrow = 100)
    rownames(m) <- paste0("G", 1:100)
    colnames(m) <- paste0("C", 1:500)
    meta <- data.frame(
        cluster = sample(c("a", "b"), 500, replace = TRUE),
        donor = sample(paste0("d", 1:6), 500, replace = TRUE)
    )

    ## collapse cells into per-(cluster, donor) pseudobulks
    dc <- collapse_counts(m, meta, c("cluster", "donor"))

    ## test each cluster against the rest
    res <- pseudobulk_deseq2(
        ~cluster,
        dc$meta_data["cluster"],
        dc$counts_mat,
        verbose = FALSE,
        present_in_min_samples = 1,
        collapse_background = FALSE,
        mode = "one_vs_all"
    )
    head(res)
}
#>   group feature baseMean log2FoldChange      lfcSE     stat     pvalue     padj
#> 1     a     G92 190.4614      0.1448955 0.07400614 1.957885 0.05024350 0.979238
#> 2     a      G6 184.2070      0.1467020 0.07623659 1.924299 0.05431713 0.979238
#> 3     a     G67 187.5035      0.1377440 0.07599174 1.812618 0.06989075 0.979238
#> 4     a     G41 190.1180      0.1324874 0.07422529 1.784937 0.07427158 0.979238
#> 5     a     G23 188.1956      0.1205691 0.08044927 1.498697 0.13395223 0.979238
#> 6     a      G8 187.4472      0.1188057 0.07962251 1.492112 0.13566991 0.979238
# }
```
