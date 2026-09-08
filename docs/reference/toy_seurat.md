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
```
