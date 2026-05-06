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
if (FALSE) { # \dontrun{
 ## generate a tiny toy dataset
 set.seed(42)
 exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
                 dimnames = list(paste0("G", 1:25), NULL))
 y <- rep(c("A", "B", "C"), each = 50)

 ## on a dense matrix
 head(wilcoxauc(exprs, y))

 ## with only some groups
 head(wilcoxauc(exprs, y, c('A', 'B')))

 ## on a sparse matrix
 exprs_sparse <- as(exprs, 'dgCMatrix')
 head(wilcoxauc(exprs_sparse, y))

 ## on a Seurat V3 object
 if (requireNamespace("Seurat", quietly = TRUE)) {
     pkg_version <- packageVersion('Seurat')
     if (pkg_version >= "3.0" & pkg_version < "4.0") {
         data(object_seurat)
         head(wilcoxauc(object_seurat, 'cell_type'))
     }
 }

 ## on a SingleCellExperiment object
 if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
     data(object_sce)
     head(wilcoxauc(object_sce, 'cell_type'))
 }
} # }
```
