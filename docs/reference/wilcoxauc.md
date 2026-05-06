# Fast Wilcoxon rank sum test and auROC

Computes auROC and Wilcoxon p-value based on Gaussian approximation.
Inputs can be

- Dense matrix or data.frame

- Sparse matrix, such as dgCMatrix

- Seurat V3 object

- SingleCellExperiment object

For detailed examples, consult the presto vignette.

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

  A feature-by-sample matrix, Seurat object, or SingleCellExperiment
  object

- ...:

  input specific parameters.

- group_by:

  (Seurat & SCE) name of groups variable ('e.g. Cluster').

- assay:

  (Seurat & SCE) name of feature matrix slot (e.g. 'data' or
  'logcounts').

- groups_use:

  (optional) which groups from y vector to test.

- seurat_assay:

  (Seurat) name of Seurat Assay (e.g. 'RNA').

- y:

  vector of group labels.

- verbose:

  boolean, TRUE for warnings and messages.

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

## Examples

``` r
if (FALSE) { # \dontrun{
 data(exprs)
 data(y)

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
