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
    ## 40 genes x 300 cells from 2 clusters across 6 donors
    m <- matrix(sample.int(8, 40 * 300, replace = TRUE), nrow = 40)
    rownames(m) <- paste0("G", 1:40)
    colnames(m) <- paste0("C", 1:300)
    meta <- data.frame(
        cluster = sample(c("a", "b"), 300, replace = TRUE),
        donor = sample(paste0("d", 1:6), 300, replace = TRUE)
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
#> 1     a      G3 113.1451      0.2054373 0.08921427 2.302740 0.02129346 0.572847
#> 2     a     G12 108.6964      0.2280651 0.10421667 2.188374 0.02864235 0.572847
#> 3     a     G13 103.3259      0.1215659 0.08240200 1.475279 0.14013755 0.982503
#> 4     a      G6 110.3473      0.1306091 0.09479678 1.377780 0.16827136 0.982503
#> 5     a      G4 106.6365      0.1401954 0.10660014 1.315152 0.18845882 0.982503
#> 6     a      G7 111.2432      0.1016056 0.09145469 1.110994 0.26657085 0.982503
# }
```
