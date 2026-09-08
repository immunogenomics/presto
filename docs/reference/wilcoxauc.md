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
wilcoxauc(
  X,
  y,
  groups_use = NULL,
  verbose = TRUE,
  nthreads = 1,
  transposed = FALSE,
  ...
)
```

## Arguments

- X:

  Input data. One of:

  - a numeric feature-by-observation matrix or `data.frame`,

  - a sparse `dgCMatrix` of the same shape,

  - a disk-backed `DelayedMatrix` (e.g. HDF5-backed, from the
    DelayedArray / HDF5Array packages), which is processed in feature
    blocks so the whole matrix never has to be loaded in memory,

  - any other matrix-like class with an `as(., "dgCMatrix")` coercion
    method (e.g. BPCells), which is converted up front,

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

- transposed:

  Set to `TRUE` when your observations (cells, samples) are in the
  **rows** of `X` and the features in the columns – i.e. `X` is the
  transpose of the default features-by-observations layout. The test
  then runs directly on that layout without materializing a transposed
  copy, which saves time and memory on large matrices. Same convention
  as the `transposed` argument of scater's `calculatePCA()` and
  `calculateUMAP()`. Only applies to matrix-like input (the `Seurat` /
  `SingleCellExperiment` dispatchers always extract
  features-by-observations). Default `FALSE`.

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
#>   feature  group  avgExpr     logFC statistic        auc         pval
#> 1      G1 jurkat 7.044320  2.402634   22061.0 0.98048889 4.798334e-47
#> 2      G2 jurkat 7.023134  2.436198   22127.0 0.98342222 1.345320e-47
#> 3      G3 jurkat 7.057136  2.092308   21980.0 0.97688889 2.524538e-46
#> 4      G4 jurkat 3.968689 -3.058986     411.0 0.01826667 1.851896e-47
#> 5      G5 jurkat 4.521172 -2.516202     496.0 0.02204444 1.392045e-46
#> 6      G6 jurkat 4.497776 -2.538578     566.5 0.02517778 5.193532e-46
#>           padj    pct_in   pct_out
#> 1 3.198889e-46 100.00000  75.33333
#> 2 1.851896e-46 100.00000  75.33333
#> 3 1.009815e-45 100.00000  80.66667
#> 4 1.851896e-46  64.00000 100.00000
#> 5 6.960227e-46  74.00000 100.00000
#> 6 1.731177e-45  72.66667 100.00000

## on a SingleCellExperiment object
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    object_sce <- toy_sce()
    head(wilcoxauc(object_sce, 'cell_type'))
}
#>   feature  group   avgExpr      logFC statistic        auc         pval
#> 1      G1 jurkat 1.6979217  0.8980784   22061.0 0.98048889 4.800228e-47
#> 2      G2 jurkat 1.6801855  0.9135012   22129.0 0.98351111 1.294566e-47
#> 3      G3 jurkat 1.7081266  0.8529468   21980.0 0.97688889 2.524883e-46
#> 4      G4 jurkat 0.6963414 -0.9875157     411.0 0.01826667 1.852288e-47
#> 5      G5 jurkat 0.7637258 -0.9285427     496.0 0.02204444 1.393166e-46
#> 6      G6 jurkat 0.7857732 -0.9053388     566.5 0.02517778 5.194357e-46
#>           padj    pct_in   pct_out
#> 1 3.200152e-46 100.00000  75.33333
#> 2 1.852288e-46 100.00000  75.33333
#> 3 1.009953e-45 100.00000  80.66667
#> 4 1.852288e-46  64.00000 100.00000
#> 5 6.965829e-46  74.00000 100.00000
#> 6 1.731452e-45  72.66667 100.00000
```
