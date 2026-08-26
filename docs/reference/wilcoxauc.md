# Fast Wilcoxon rank-sum test and auROC across groups

For every (feature, group) pair, computes the Wilcoxon rank-sum
statistic comparing observations in that group against all other
observations, and the area under the ROC curve as a measure of
separability. P-values come from the standard Gaussian approximation to
the U statistic with a tie correction. Returns one row per (feature,
group) with effect-size and percent-expressed columns alongside the test
statistics.

## Usage

``` r
wilcoxauc(X, ...)

# S3 method for class 'seurat'
wilcoxauc(X, ...)

# S3 method for class 'Seurat'
wilcoxauc(
  X,
  group_by = NULL,
  assay = "data",
  groups_use = NULL,
  seurat_assay = "RNA",
  ...
)

# S3 method for class 'SingleCellExperiment'
wilcoxauc(X, group_by = NULL, assay = NULL, groups_use = NULL, ...)

# Default S3 method
wilcoxauc(X, y, groups_use = NULL, verbose = TRUE, nthreads = 1, ...)
```

## Arguments

- X:

  Input data. One of:

  - a numeric feature-by-observation matrix or `data.frame`,

  - a sparse `dgCMatrix` of the same shape,

  - a `Seurat` (v3+) object,

  - a `SingleCellExperiment` object.

  `X` must not contain `NA` values. Unlike
  [`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html),
  which drops missing values per observation, `wilcoxauc()` errors on
  `NA` input rather than returning silently incorrect results; remove or
  impute missing values first.

- ...:

  Passed to the input-specific method.

- group_by:

  For `Seurat` and `SingleCellExperiment` input, name of the metadata
  column that holds the group labels (e.g. `"cluster"`). For `Seurat`,
  defaults to `Idents(X)`.

- assay:

  For `Seurat`, the layer name within the selected assay (e.g. `"data"`,
  `"counts"`, `"scale.data"`). For `SingleCellExperiment`, the assay
  name (e.g. `"logcounts"`, `"counts"`). Defaults pick a sensible value
  per input class.

- groups_use:

  Optional character vector restricting the test to a subset of groups
  in `y` (or `group_by`). Default `NULL` tests every group.

- seurat_assay:

  For `Seurat` input, the name of the assay to pull from (e.g. `"RNA"`).
  Default `"RNA"`.

- y:

  For matrix input, a character/factor vector of group labels with
  length equal to `ncol(X)`. Ignored for `Seurat` /
  `SingleCellExperiment` input (use `group_by` instead).

- verbose:

  Logical. Print warnings and informational messages. Default `TRUE`.

- nthreads:

  Number of threads for the per-feature ranking of sparse (`dgCMatrix`)
  input. Default `1` (serial). Values above `1` split the ranking across
  threads; the result is identical regardless of the thread count. Only
  sparse input is parallelized – dense matrix and `data.frame` input are
  always processed serially. When running under `R CMD check` or on
  CRAN, keep this at the default so no more than two cores are used.

## Value

table with the following columns:

- **feature** - feature name (e.g. gene name).

- **group** - group name.

- **avgExpr** - mean value of feature in group.

- **logFC** - difference of mean feature values between observations in
  the group vs out of the group. When the input is log-transformed
  expression (e.g. Seurat's `"data"` layer or `logcounts`), this
  difference of means is a log fold change. On raw (untransformed)
  values it is a plain difference of means, not a fold change.

- **statistic** - Wilcoxon rank sum U statistic.

- **auc** - area under the receiver operator curve.

- **pval** - nominal p value.

- **padj** - Benjamini-Hochberg adjusted p value.

- **pct_in** - Percent of observations in the group with non-zero
  feature value.

- **pct_out** - Percent of observations out of the group with non-zero
  feature value.

## Details

Designed to be fast enough to run on whole-genome × hundred-thousand-
cell single-cell matrices in seconds. Sparse `dgCMatrix` inputs are
processed without densification. Convenience dispatchers extract the
counts matrix and group labels from `Seurat` and `SingleCellExperiment`
objects. See the `getting-started` vignette for an end-to-end example on
a real dataset.

## See also

[`top_markers()`](https://immunogenomics.github.io/presto/reference/top_markers.md)
to summarize markers per group;
[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md)
for a count-based pseudobulk alternative.

## Examples

``` r
## generate a tiny toy dataset
set.seed(42)
exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
                dimnames = list(paste0("G", 1:25), NULL))
y <- rep(c("A", "B", "C"), each = 50)

## on a dense matrix
head(wilcoxauc(exprs, y))
#>   feature group avgExpr logFC statistic    auc       pval      padj pct_in
#> 1      G1     A    2.10  0.20    2740.0 0.5480 0.32643797 0.5100593     86
#> 2      G2     A    1.58 -0.56    1896.0 0.3792 0.01372592 0.1715740     84
#> 3      G3     A    1.86  0.02    2432.5 0.4865 0.78263792 0.8933109     84
#> 4      G4     A    1.96 -0.21    2384.0 0.4768 0.63616986 0.8835693     90
#> 5      G5     A    2.00 -0.28    2214.0 0.4428 0.24419173 0.5052190     82
#> 6      G6     A    2.26  0.24    2775.0 0.5550 0.26271390 0.5052190     90
#>   pct_out
#> 1      89
#> 2      86
#> 3      89
#> 4      90
#> 5      88
#> 6      83

## restrict the comparison to a subset of groups
head(wilcoxauc(exprs, y, c('A', 'B')))
#>   feature group avgExpr logFC statistic    auc       pval      padj pct_in
#> 1      G1     A    2.10  0.02      1284 0.5136 0.81297801 0.9680542     86
#> 2      G2     A    1.58 -0.52       950 0.3800 0.03365247 0.4716572     84
#> 3      G3     A    1.86 -0.02      1197 0.4788 0.70869519 0.9324937     84
#> 4      G4     A    1.96 -0.28      1158 0.4632 0.51725570 0.8474382     90
#> 5      G5     A    2.00 -0.02      1236 0.4944 0.92392119 0.9802095     82
#> 6      G6     A    2.26  0.36      1437 0.5748 0.18870167 0.6739345     90
#>   pct_out
#> 1      90
#> 2      86
#> 3      90
#> 4      88
#> 5      84
#> 6      82

## on a sparse matrix
exprs_sparse <- as(exprs, 'dgCMatrix')
head(wilcoxauc(exprs_sparse, y))
#>   feature group avgExpr logFC statistic    auc       pval      padj pct_in
#> 1      G1     A    2.10  0.20    2740.0 0.5480 0.32643797 0.5100593     86
#> 2      G2     A    1.58 -0.56    1896.0 0.3792 0.01372592 0.1715740     84
#> 3      G3     A    1.86  0.02    2432.5 0.4865 0.78263792 0.8933109     84
#> 4      G4     A    1.96 -0.21    2384.0 0.4768 0.63616986 0.8835693     90
#> 5      G5     A    2.00 -0.28    2214.0 0.4428 0.24419173 0.5052190     82
#> 6      G6     A    2.26  0.24    2775.0 0.5550 0.26271390 0.5052190     90
#>   pct_out
#> 1      89
#> 2      86
#> 3      89
#> 4      90
#> 5      88
#> 6      83

## on a Seurat object (>= v3)
if (requireNamespace("Seurat", quietly = TRUE) &&
    packageVersion("Seurat") >= "3.0") {
    object_seurat <- toy_seurat()
    head(wilcoxauc(object_seurat, 'cell_type'))
}
#>   feature  group  avgExpr     logFC statistic       auc         pval
#> 1      G1 jurkat 6.827790  1.657491   20418.5 0.9074889 2.870368e-34
#> 2      G2 jurkat 6.816695  1.517735   19311.5 0.8582889 7.178349e-27
#> 3      G3 jurkat 6.781403  1.482781   20330.5 0.9035778 1.213350e-33
#> 4      G4 jurkat 5.388072 -1.357031    3071.0 0.1364889 1.306183e-27
#> 5      G5 jurkat 5.150948 -1.581087    2800.0 0.1244444 2.298078e-29
#> 6      G6 jurkat 5.282571 -1.406932    3370.0 0.1497778 9.540231e-26
#>           padj    pct_in  pct_out
#> 1 5.740736e-33 100.00000 86.00000
#> 2 2.871340e-26 100.00000 86.66667
#> 3 1.213350e-32  98.66667 88.66667
#> 4 6.530914e-27  88.66667 98.66667
#> 5 1.532052e-28  85.33333 99.33333
#> 6 3.180077e-25  88.00000 98.66667

## on a SingleCellExperiment object
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    object_sce <- toy_sce()
    head(wilcoxauc(object_sce, 'cell_type'))
}
#>   feature  group   avgExpr      logFC statistic       auc         pval
#> 1      G1 jurkat 1.7244684  0.7575272   20420.5 0.9075778 2.778842e-34
#> 2      G2 jurkat 1.7182759  0.6798860   19310.5 0.8582444 7.284606e-27
#> 3      G3 jurkat 1.7384917  0.7613587   20330.0 0.9035556 1.223619e-33
#> 4      G4 jurkat 1.0381649 -0.6710793    3070.5 0.1364667 1.296878e-27
#> 5      G5 jurkat 0.9718625 -0.7058267    2798.5 0.1243778 2.246744e-29
#> 6      G6 jurkat 0.9863318 -0.6863914    3369.5 0.1497556 9.473597e-26
#>           padj    pct_in  pct_out
#> 1 5.557684e-33 100.00000 86.00000
#> 2 2.913842e-26 100.00000 86.66667
#> 3 1.223619e-32  98.66667 88.66667
#> 4 6.484390e-27  88.66667 98.66667
#> 5 1.497830e-28  85.33333 99.33333
#> 6 3.157866e-25  88.00000 98.66667
```
