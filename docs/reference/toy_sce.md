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
