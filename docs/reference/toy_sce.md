# Toy SingleCellExperiment object for examples and tests

Builds a small deterministic `SingleCellExperiment` (20 genes x 300
cells, two cell types in the `cell_type` colData column) with the
installed SingleCellExperiment version, so no serialized object needs to
ship with the package. The same simulated counts as
[`toy_seurat()`](https://immunogenomics.github.io/presto/reference/toy_seurat.md)
are stored in the `counts` assay, with log-normalized values in
`logcounts`. The generator restores the caller's random-number state.

## Usage

``` r
toy_sce()
```

## Value

A `SingleCellExperiment` with `counts` and `logcounts` assays and a
`cell_type` colData column with values `"jurkat"` and `"t293"`.

## Details

Requires the suggested package `SingleCellExperiment` (Bioconductor).

## See also

[`toy_seurat()`](https://immunogenomics.github.io/presto/reference/toy_seurat.md),
[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)

## Examples

``` r
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    object_sce <- toy_sce()
    head(wilcoxauc(object_sce, "cell_type"))
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
