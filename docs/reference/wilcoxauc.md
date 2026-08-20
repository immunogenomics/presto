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
wilcoxauc(X, y, groups_use = NULL, verbose = TRUE, ...)
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

## Value

table with the following columns:

- **feature** - feature name (e.g. gene name).

- **group** - group name.

- **avgExpr** - mean value of feature in group.

- **logFC** - log fold change between observations in group vs out.

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
#> 1      G1     A    2.10  0.20    2740.0 0.5480 0.32643843 0.5100600     86
#> 2      G2     A    1.58 -0.56    1896.0 0.3792 0.01372592 0.1715740     84
#> 3      G3     A    1.86  0.02    2432.5 0.4865 0.78264917 0.8933112     84
#> 4      G4     A    1.96 -0.21    2384.0 0.4768 0.63616986 0.8835693     90
#> 5      G5     A    2.00 -0.28    2214.0 0.4428 0.24419217 0.5052190     82
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
#> 1      G1     A    2.10  0.02      1284 0.5136 0.81297858 0.9680542     86
#> 2      G2     A    1.58 -0.52       950 0.3800 0.03365247 0.4716572     84
#> 3      G3     A    1.86 -0.02      1197 0.4788 0.70871296 0.9325171     84
#> 4      G4     A    1.96 -0.28      1158 0.4632 0.51725570 0.8474382     90
#> 5      G5     A    2.00 -0.02      1236 0.4944 0.92392119 0.9802096     82
#> 6      G6     A    2.26  0.36      1437 0.5748 0.18871553 0.6739840     90
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
#> 1      G1     A    2.10  0.20    2740.0 0.5480 0.32643843 0.5100600     86
#> 2      G2     A    1.58 -0.56    1896.0 0.3792 0.01372592 0.1715740     84
#> 3      G3     A    1.86  0.02    2432.5 0.4865 0.78264917 0.8933112     84
#> 4      G4     A    1.96 -0.21    2384.0 0.4768 0.63616986 0.8835693     90
#> 5      G5     A    2.00 -0.28    2214.0 0.4428 0.24419217 0.5052190     82
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
    data(object_seurat)
    object_seurat <- Seurat::UpdateSeuratObject(object_seurat)
    head(wilcoxauc(object_seurat, 'cell_type'))
}
#> Validating object structure
#> Updating object slots
#> Ensuring keys are in the proper structure
#> Updating matrix keys for DimReduc ‘pca’
#> Ensuring keys are in the proper structure
#> Ensuring feature names don't have underscores or pipes
#> Updating slots in RNA
#> Updating slots in pca
#> Setting assay used for NormalizeData.RNA to RNA
#> Setting assay used for FindVariableFeatures.RNA to RNA
#> Setting assay used for ScaleData.RNA to RNA
#> Setting assay used for RunPCA.RNA to RNA
#> Validating object structure for Assay ‘RNA’
#> Validating object structure for DimReduc ‘pca’
#> Object representation is consistent with the most current Seurat version
#>   feature  group  avgExpr        logFC statistic       auc      pval      padj
#> 1      G1 jurkat 5.957396  0.040642962     11343 0.5043351 0.8972436 0.9989377
#> 2      G2 jurkat 5.945143 -0.027840659     11209 0.4983771 0.9617722 0.9989377
#> 3      G3 jurkat 5.919214 -0.058395057     10306 0.4582277 0.2112380 0.9913131
#> 4      G4 jurkat 5.931007 -0.059677089     10859 0.4828153 0.6073122 0.9989377
#> 5      G5 jurkat 5.957468 -0.020184665     11244 0.4999333 0.9989377 0.9989377
#> 6      G6 jurkat 5.919087 -0.007454977     11038 0.4907741 0.7828582 0.9989377
#>   pct_in pct_out
#> 1    100     100
#> 2    100     100
#> 3    100     100
#> 4    100     100
#> 5    100     100
#> 6    100     100

## on a SingleCellExperiment object
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    data(object_sce)
    head(wilcoxauc(object_sce, 'cell_type'))
}
#>   feature  group  avgExpr        logFC statistic       auc      pval      padj
#> 1      G1 jurkat 5.957396  0.040642962     11343 0.5043351 0.8972436 0.9989377
#> 2      G2 jurkat 5.945143 -0.027840659     11209 0.4983771 0.9617722 0.9989377
#> 3      G3 jurkat 5.919214 -0.058395057     10306 0.4582277 0.2112380 0.9913131
#> 4      G4 jurkat 5.931007 -0.059677089     10859 0.4828153 0.6073122 0.9989377
#> 5      G5 jurkat 5.957468 -0.020184665     11244 0.4999333 0.9989377 0.9989377
#> 6      G6 jurkat 5.919087 -0.007454977     11038 0.4907741 0.7828582 0.9989377
#>   pct_in pct_out
#> 1    100     100
#> 2    100     100
#> 3    100     100
#> 4    100     100
#> 5    100     100
#> 6    100     100
```
