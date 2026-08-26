# Changelog

## presto 1.1.0

First CRAN release. Version 1.0.0 refers to the pre-CRAN package
installed from GitHub.

### New features

- [`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
  gains an `nthreads` argument: the per-feature ranking of sparse
  (`dgCMatrix`) input can run multithreaded. Results are identical at
  any thread count; roughly 2.6x faster with 4 threads on a 3,000 x
  80,000 matrix. The default (`nthreads = 1`) is unchanged.
- New
  [`toy_seurat()`](https://immunogenomics.github.io/presto/reference/toy_seurat.md)
  and
  [`toy_sce()`](https://immunogenomics.github.io/presto/reference/toy_sce.md)
  generate small example Seurat / SingleCellExperiment objects on the
  fly with the installed package versions, so they can never go stale.
- New
  [`load_ircolitis_cd8()`](https://immunogenomics.github.io/presto/reference/load_ircolitis_cd8.md)
  downloads and caches the 25,341-cell colon CD8 T-cell dataset
  (GSE206299, Thomas et al. 2024) used by the vignettes.
- Two vignettes: *Getting started with presto* and *Pseudobulk
  differential expression with DESeq2*, both run on the real ircolitis
  dataset.

### Breaking changes

- [`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
  now stops with an informative error when `X` contains `NA` values,
  instead of silently returning incorrect results (#25).
- The shipped `object_seurat` and `object_sce` datasets were removed;
  use
  [`toy_seurat()`](https://immunogenomics.github.io/presto/reference/toy_seurat.md)
  and
  [`toy_sce()`](https://immunogenomics.github.io/presto/reference/toy_sce.md)
  instead.
- `Rcpp` and `data.table` moved from `Depends` to `Imports`, so they are
  no longer attached to the search path when presto is loaded.

### Bug fixes

- The Wilcoxon variance correction now includes every tie group: the
  final tie group per feature and the implicit-zero group of all-zero
  features used to be dropped, inflating p-values on heavily tied data
  (#29). P-values now agree with
  [`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) to
  machine precision; on typical count data the change is tiny (\< 1e-4),
  but it is meaningful for binary or heavily tied features. Constant
  features report `pval = 1` instead of an arbitrary value.
- [`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
  and
  [`rank_matrix()`](https://immunogenomics.github.io/presto/reference/rank_matrix.md)
  no longer modify their input in place: dense input used to be silently
  overwritten with ranks via an aliased memory buffer (#7). Dense
  `avgExpr` and `logFC` are computed from the original values (also
  \#7).
- The Seurat dispatcher of
  [`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
  uses the `layer` argument required by Seurat 5 (#44).
- [`pseudobulk_within()`](https://immunogenomics.github.io/presto/reference/pseudobulk_within.md)
  handles character and multi-level contrast variables.
- [`pseudobulk_pairwise()`](https://immunogenomics.github.io/presto/reference/pseudobulk_pairwise.md)
  no longer fails on single-column `meta_data`.
- [`summarize_dge_pairs()`](https://immunogenomics.github.io/presto/reference/summarize_dge_pairs.md)
  no longer prints debug output, and
  [`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md)
  no longer warns unconditionally.

### Performance

- The sparse-input statistics of
  [`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
  (transpose, ranking, and all group-wise reductions) are fused into a
  single-pass C++ kernel, about 12% faster serially before any
  threading.

## presto 1.0.0

- Original GitHub release.
