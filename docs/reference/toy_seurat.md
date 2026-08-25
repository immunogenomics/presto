# Toy Seurat object for examples and tests

Builds a small deterministic `Seurat` object (20 genes x 300 cells, two
cell types in the `cell_type` metadata column) with the installed Seurat
version, so no serialized object needs to ship with the package or be
updated when Seurat changes its internal class structure. Counts are
Poisson-simulated with a few upregulated marker genes per cell type, and
log-normalized into the `data` layer. The generator restores the
caller's random-number state, so calling it does not perturb
reproducibility.

## Usage

``` r
toy_seurat()
```

## Value

A `Seurat` object with `counts` and `data` layers and a `cell_type`
metadata column with values `"jurkat"` and `"t293"`.

## Details

Requires the suggested package `Seurat`.

## See also

[`toy_sce()`](https://immunogenomics.github.io/presto/reference/toy_sce.md),
[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE)) {
    object_seurat <- toy_seurat()
    head(wilcoxauc(object_seurat, "cell_type"))
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
```
