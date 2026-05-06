# presto

Fast differential expression and marker discovery for single-cell data.

`presto` is built around two function families:

- **`wilcoxauc()`** — Wilcoxon rank-sum test and area-under-the-ROC
  across groups, fast enough to run on whole-genome × hundred-thousand-
  cell matrices in seconds.
- **`pseudobulk_deseq2()`** — collapse cells into per-donor pseudobulks
  with `collapse_counts()` and run a count-based regression with
  `DESeq2` in one-vs-all, pairwise, or within-cluster designs.

Both work directly with sparse `dgCMatrix` input and have dispatchers
for [`Seurat`](https://satijalab.org/seurat/) and
[`SingleCellExperiment`](https://bioconductor.org/packages/SingleCellExperiment/)
objects.

## Installation

```r
# install.packages("devtools")
devtools::install_github("immunogenomics/presto")
```

`pseudobulk_deseq2()` additionally requires
[`DESeq2`](https://bioconductor.org/packages/DESeq2/), and the
`load_ircolitis_cd8()` demo loader requires `rhdf5`. Both live in
`Suggests`, so you only need them if you use those features.

## Quick example

```r
library(presto)

data(exprs)  # 25-gene × 150-cell toy matrix
data(y)      # 3 group labels

res <- wilcoxauc(exprs, y)
head(res)
#>   feature group avgExpr   logFC statistic   auc   pval  padj pct_in pct_out
#> 1      G1     A   1.145  0.0823      2890 0.553 0.2572 0.643   74.5    58.9
#> 2      G2     A   0.764 -0.2995      2190 0.419 0.0809 0.409   49.1    65.3
#> 3      G3     A   0.982 -0.1866      2429 0.465 0.4522 0.748   65.5    68.4
#> 4      G4     A   1.218  0.3971      3176 0.608 0.0200 0.409   70.9    55.8
#> 5      G5     A   1.418  0.3656      3026 0.579 0.0931 0.409   72.7    56.8
#> 6      G6     A   0.927 -0.1780      2424 0.464 0.4387 0.748   56.4    60.0
```

`top_markers()` summarises the most distinguishing features per group:

```r
top_markers(res, n = 5, auc_min = 0.6)
#> # A tibble: 5 × 4
#>    rank A     B     C
#>   <int> <chr> <chr> <chr>
#> 1     1 G4    G20   G1
#> 2     2 NA    G5    G21
#> 3     3 NA    G15   G6
#> 4     4 NA    G19   G16
#> 5     5 NA    G25   G11
```

The same call also works on Seurat and SingleCellExperiment objects:

```r
wilcoxauc(seurat_object, group_by = "cluster")
wilcoxauc(sce_object,    group_by = "cluster")
```

## Walkthroughs

For full examples on a real 25,000-cell single-cell RNA-seq dataset, see:

- [**Getting started**](https://immunogenomics.github.io/presto/articles/getting-started.html)
  — `wilcoxauc()` and `top_markers()` for marker discovery, including a
  description of every output column and how to restrict the comparison
  to a subset of groups.
- [**Pseudobulk differential expression with DESeq2**](https://immunogenomics.github.io/presto/articles/pseudobulk.html)
  — `collapse_counts()` and the three `pseudobulk_deseq2()` modes
  (one-vs-all, pairwise, within), plus `top_markers_dds()` and
  `summarize_dge_pairs()`.

## Performance

A benchmark on 1,000,000 observations × 1,000 features × 10 groups
runs in 16 seconds on sparse input and 85 seconds on dense input.
