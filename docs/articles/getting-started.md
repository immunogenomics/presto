# Getting started with presto

`presto` runs the Wilcoxon rank-sum test and computes the auROC
statistic quickly enough to be used as a routine marker-finding step on
whole-genome, hundred-thousand-cell single-cell experiments.

## Demo data

This vignette uses CD8 T cells from human colon tissue published by
[Thomas et al., *Nature Medicine*
2024](https://www.nature.com/articles/s41591-024-02895-x)
([GSE206299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206299)).
[`load_ircolitis_cd8()`](https://immunogenomics.github.io/presto/reference/load_ircolitis_cd8.md)
downloads the 121 MB `.h5ad.gz` file on first use, caches it under
`tools::R_user_dir("presto", "cache")`, and parses it with `rhdf5`.

``` r

library(presto)
library(Matrix)
library(dplyr)

d <- load_ircolitis_cd8(verbose = FALSE)
dim(d$counts)
#> [1] 28165 25341
```

`d$counts` is a 28,165 × 25,341 sparse matrix of UMI counts (genes ×
cells). `d$obs` carries cell metadata, including a `cluster` column we
will use to identify markers.

``` r

table(d$obs$cluster)
#> 
#>    1    2    3    4    5    6    7    8    9 
#> 8380 5115 2604 2246 2235 2049 1860  655  197
```

## Quick start

The simplest call passes a feature-by-cell matrix and a vector of group
labels:

``` r

res <- wilcoxauc(d$counts, d$obs$cluster)
head(res)
#>                       feature group  avgExpr     logFC statistic   auc     pval    padj pct_in
#> 1 ENSG00000243485|MIR1302-2HG     1 0.000000  0.000000  71066590 0.500 1.000000 1.00000 0.0000
#> 2  ENSG00000238009|AL627309.1     1 0.000239 -0.000233  71050031 0.500 0.379630 0.67361 0.0239
#> 3  ENSG00000241599|AL627309.4     1 0.000000  0.000000  71066590 0.500 1.000000 1.00000 0.0000
#> 4  ENSG00000235146|AC114498.1     1 0.000000  0.000000  71066590 0.500 1.000000 1.00000 0.0000
#> 5  ENSG00000237491|AL669831.5     1 0.053938  0.011016  71757724 0.505 0.000418 0.00187 5.1074
#> 6      ENSG00000177757|FAM87B     1 0.000000  0.000000  71066590 0.500 1.000000 1.00000 0.0000
#>   pct_out
#> 1  0.0000
#> 2  0.0472
#> 3  0.0000
#> 4  0.0000
#> 5  4.1389
#> 6  0.0000
```

Each row reports one (gene, group) test. Rows are returned for every
gene × cluster combination, comparing observations in that cluster
against all other observations.

## Result columns

| column      | description                                                   |
|-------------|---------------------------------------------------------------|
| `feature`   | name of the feature (gene)                                    |
| `group`     | name of the group / cluster                                   |
| `avgExpr`   | mean expression of the feature in the group                   |
| `logFC`     | log fold change between in-group vs out-of-group observations |
| `statistic` | Wilcoxon rank-sum U statistic                                 |
| `auc`       | area under the ROC curve                                      |
| `pval`      | nominal p-value, two-sided Gaussian approximation of U        |
| `padj`      | Benjamini-Hochberg adjusted p-value                           |
| `pct_in`    | percent of in-group cells with non-zero expression            |
| `pct_out`   | percent of out-of-group cells with non-zero expression        |

## Top markers per cluster

[`top_markers()`](https://immunogenomics.github.io/presto/reference/top_markers.md)
extracts the most distinguishing features per group. Use `auc_min`,
`padj_max`, `pct_in_min`, and `pct_out_max` to filter.

``` r

top_markers(res, n = 5, auc_min = 0.6, padj_max = 1e-3)
#> # A tibble: 5 × 10
#>    rank `1`                    `2`                      `3`      `4`   `5`   `6`   `7`   `8`   `9`  
#>   <int> <chr>                  <chr>                    <chr>    <chr> <chr> <chr> <chr> <chr> <chr>
#> 1     1 ENSG00000198727|MT-CYB ENSG00000153563|CD8A     ENSG000… ENSG… ENSG… ENSG… ENSG… ENSG… ENSG…
#> 2     2 ENSG00000117632|STMN1  ENSG00000142669|SH3BGRL3 ENSG000… ENSG… ENSG… ENSG… ENSG… ENSG… ENSG…
#> 3     3 ENSG00000128951|DUT    ENSG00000196154|S100A4   ENSG000… ENSG… ENSG… ENSG… ENSG… ENSG… ENSG…
#> 4     4 ENSG00000198804|MT-CO1 ENSG00000116824|CD2      ENSG000… ENSG… ENSG… ENSG… ENSG… ENSG… ENSG…
#> 5     5 ENSG00000198938|MT-CO3 ENSG00000167286|CD3D     ENSG000… ENSG… ENSG… ENSG… ENSG… ENSG… ENSG…
```

The gene identifiers are the verbatim row names from the input matrix
(`ENSG…|SYMBOL`). To plot just the symbol portion, split at `|`:

``` r

top5 <- top_markers(res, n = 5, auc_min = 0.6, padj_max = 1e-3)
top5 %>%
    dplyr::mutate(across(-rank, ~ sub(".*\\|", "", .x)))
#> # A tibble: 5 × 10
#>    rank `1`    `2`      `3`   `4`   `5`    `6`      `7`      `8`     `9`   
#>   <int> <chr>  <chr>    <chr> <chr> <chr>  <chr>    <chr>    <chr>   <chr> 
#> 1     1 MT-CYB CD8A     GZMK  GAPDH IL7R   ACTG1    TYROBP   TNFRSF4 JCHAIN
#> 2     2 STMN1  SH3BGRL3 CD27  GZMA  RPS12  CCL5     FCER1G   SRGN    IGHA1 
#> 3     3 DUT    S100A4   CST7  CD7   S100A4 TMSB4X   TRDC     ITGB2   IGKC  
#> 4     4 MT-CO1 CD2      ITGB2 TPI1  LTB    CD52     TNFRSF18 NFKBIA  IGLC2 
#> 5     5 MT-CO3 CD3D     GZMH  GZMB  RPS14  SH3BGRL3 GSTP1    CTLA4   IGLC3
```

## Restricting to a subset of groups

Use `groups_use` to compare only a subset of groups against each other:

``` r

res_12 <- wilcoxauc(d$counts, d$obs$cluster, groups_use = c("1", "2"))
top_markers(res_12, n = 10, auc_min = 0.6, padj_max = 1e-3)
#> # A tibble: 10 × 3
#>     rank `1`                     `2`                     
#>    <int> <chr>                   <chr>                   
#>  1     1 ENSG00000198727|MT-CYB  ENSG00000153563|CD8A    
#>  2     2 ENSG00000198938|MT-CO3  ENSG00000116824|CD2     
#>  3     3 ENSG00000198886|MT-ND4  ENSG00000196154|S100A4  
#>  4     4 ENSG00000198712|MT-CO2  ENSG00000111796|KLRB1   
#>  5     5 ENSG00000198804|MT-CO1  ENSG00000197956|S100A6  
#>  6     6 ENSG00000117632|STMN1   ENSG00000142669|SH3BGRL3
#>  7     7 ENSG00000212907|MT-ND4L ENSG00000137078|SIT1    
#>  8     8 ENSG00000198888|MT-ND1  ENSG00000166710|B2M     
#>  9     9 ENSG00000128951|DUT     ENSG00000213658|LAT     
#> 10    10 ENSG00000164104|HMGB2   ENSG00000197540|GZMM
```

## Seurat and SingleCellExperiment

[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
dispatches on the input class. For Seurat or SingleCellExperiment
objects, pass the object and the name of the metadata column that
defines the groups:

``` r

# Seurat
wilcoxauc(seurat_object, group_by = "cluster")

# SingleCellExperiment
wilcoxauc(sce_object, group_by = "cluster")
```

For the full set of arguments (assay selection, layer, etc.), see
[`?wilcoxauc`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md).

## Notes on input types

[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
accepts dense matrices, sparse `dgCMatrix` inputs, and `data.frame`s.
Sparse input is faster and uses far less memory on single-cell data. The
matrix should be **features × observations** (rows are genes, columns
are cells).

## Session info

``` r

sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Sequoia 15.6.1
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US/en_US/en_US/C/en_US/en_US
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dplyr_1.2.1   Matrix_1.7-4  presto_1.0.0  knitr_1.51    ggplot2_4.0.2
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.2         cli_3.6.5           rlang_1.1.7         xfun_0.56          
#>  [5] otel_0.2.0          purrr_1.2.1         generics_0.1.4      S7_0.2.1           
#>  [9] data.table_1.18.2.1 glue_1.8.0          scales_1.4.0        grid_4.5.2         
#> [13] evaluate_1.0.5      tibble_3.3.1        Rhdf5lib_1.32.0     lifecycle_1.0.5    
#> [17] compiler_4.5.2      RColorBrewer_1.1-3  rhdf5filters_1.22.0 Rcpp_1.1.1         
#> [21] pkgconfig_2.0.3     tidyr_1.3.2         rhdf5_2.54.1        lattice_0.22-7     
#> [25] farver_2.1.2        R6_2.6.1            utf8_1.2.6          dichromat_2.0-0.1  
#> [29] tidyselect_1.2.1    pillar_1.11.1       magrittr_2.0.5      tools_4.5.2        
#> [33] withr_3.0.2         gtable_0.3.6
```
