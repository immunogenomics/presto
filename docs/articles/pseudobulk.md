# Pseudobulk differential expression with DESeq2

## Why pseudobulk?

Single cells from the same donor are not statistically independent: they
share genetic background, batch effects, and ambient RNA. Treating each
cell as an independent observation in a regression inflates the
effective sample size and distorts p-values.

A common, principled remedy is **pseudobulk**: sum each donor’s cells
within a group (e.g. cluster) to produce a per-donor count vector, then
run a count-based regression like
[`DESeq2::DESeq`](https://rdrr.io/pkg/DESeq2/man/DESeq.html) over
donors. The unit of replication becomes the donor, not the cell.

`presto` provides:

- [`collapse_counts()`](https://immunogenomics.github.io/presto/reference/collapse_counts.md)
  — fast collapse of a cell-by-gene count matrix to pseudobulks defined
  by one or more metadata columns.
- [`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md)
  — run DESeq2 across the resulting pseudobulks in a one-vs-all,
  pairwise, or within-group design.
- [`top_markers_dds()`](https://immunogenomics.github.io/presto/reference/top_markers_dds.md)
  — extract top features per group from the result.

## Demo data

We use the colon-tissue CD8 T-cell dataset from [Thomas et al. *Nat.
Med.* 2024](https://www.nature.com/articles/s41591-024-02895-x)
([GSE206299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206299)).
27 donors, 9 cell clusters, 25,341 cells.

``` r

library(presto)
library(Matrix)
library(dplyr)

d <- load_ircolitis_cd8(verbose = FALSE)
dim(d$counts)
#> [1] 28165 25341
table(d$obs$donor, d$obs$cluster)[1:5, ]
#>          
#>             1   2   3   4   5   6   7   8   9
#>   MC_1    323 110  46 100  80  51  61   7   5
#>   MC_2    580 229 113 154 179 115  78  27   5
#>   MC_9    698 449 143 174 230 119  91  32  10
#>   SIC_100 428 373 119 152 150 164  59  25  10
#>   SIC_109 172  85  33  68  45  21 170   6   4
```

The crosstab confirms most donors contribute cells to most clusters —
good news for pseudobulk, since clusters with few represented donors
will have low statistical power.

## Step 1: collapse cells to pseudobulks

[`collapse_counts()`](https://immunogenomics.github.io/presto/reference/collapse_counts.md)
sums the columns of a count matrix according to one or more metadata
columns. Here we collapse by `(donor, cluster)`, drop any (donor,
cluster) pair with fewer than 10 cells, and keep the cell count `N` so
we can inspect what came through.

``` r

md <- data.frame(
    donor   = as.character(d$obs$donor),
    cluster = as.character(d$obs$cluster),
    case    = as.character(d$obs$case),
    stringsAsFactors = FALSE
)
cc <- collapse_counts(
    d$counts, md, c("donor", "cluster", "case"),
    min_cells_per_group = 10,
    keep_n = TRUE
)
dim(cc$counts_mat)
#> [1] 28165   204
head(cc$meta_data)
#>   donor cluster    case   N
#> 1  MC_1       5 Control  80
#> 2  MC_1       3 Control  46
#> 3  MC_1       4 Control 100
#> 4  MC_1       1 Control 323
#> 5  MC_1       6 Control  51
#> 6  MC_1       2 Control 110
```

We collapse by `(donor, cluster, case)` instead of just
`(donor, cluster)` because `case` is constant within donor — including
it now makes it available for the within-cluster contrast in [Step
4](#step-4-within-cluster-case-vs-control).

`cc$counts_mat` is now a 28,165-gene × 204-pseudobulk matrix. Each
column is one (donor, cluster) sample; `cc$meta_data` is the matching
`data.frame` of pseudobulk identities.

``` r

table(cc$meta_data$cluster)
#> 
#>  1  2  3  4  5  6  7  8  9 
#> 27 27 26 27 25 23 26 18  5
```

## Step 2: one-vs-all DESeq2

[`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md)
with `mode = "one_vs_all"` tests each cluster against the union of the
others. With `collapse_background = FALSE` we keep each background
pseudobulk as its own observation; setting it to `TRUE` collapses the
background pseudobulks per donor first, which is preferable when donors
are unbalanced across clusters.

The first variable in the formula is the contrast variable. `meta_data`
should contain only the pseudobulk identifying variables used in the
formula.

``` r

md_dge <- cc$meta_data[, "cluster", drop = FALSE]
res <- pseudobulk_deseq2(
    ~cluster, md_dge, cc$counts_mat,
    verbose = FALSE,
    min_counts_per_sample = 10,
    present_in_min_samples = 5,
    collapse_background = FALSE,
    mode = "one_vs_all"
)
head(res)
#>   group                  feature baseMean log2FoldChange lfcSE stat   pvalue
#> 1     1      ENSG00000128951|DUT    42.72           1.35 0.108 12.5 6.67e-36
#> 2     1     ENSG00000105486|LIG1     4.18           1.69 0.138 12.2 2.49e-34
#> 3     1 ENSG00000134291|TMEM106C    13.14           1.67 0.138 12.1 7.54e-34
#> 4     1    ENSG00000188486|H2AFX    10.49           1.86 0.155 12.0 2.92e-33
#> 5     1     ENSG00000120802|TMPO    14.95           1.34 0.113 11.9 1.13e-32
#> 6     1     ENSG00000163535|SGO2     2.85           2.17 0.186 11.7 1.50e-31
#>       padj
#> 1 6.69e-32
#> 2 1.25e-30
#> 3 2.52e-30
#> 4 7.33e-30
#> 5 2.27e-29
#> 6 2.50e-28
```

Result columns come straight from
[`DESeq2::results()`](https://rdrr.io/pkg/DESeq2/man/results.html),
prepended with the group identifier:

| column           | description                         |
|------------------|-------------------------------------|
| `group`          | foreground cluster                  |
| `feature`        | gene identifier                     |
| `baseMean`       | mean of normalized counts           |
| `log2FoldChange` | effect size, group vs background    |
| `lfcSE`          | standard error of `log2FoldChange`  |
| `stat`           | Wald statistic                      |
| `pvalue`         | nominal p-value                     |
| `padj`           | Benjamini-Hochberg adjusted p-value |

## Step 3: top markers per cluster

[`top_markers_dds()`](https://immunogenomics.github.io/presto/reference/top_markers_dds.md)
filters by `padj_max` and `lfc_min` and returns a wide table of the top
features per group:

``` r

top5 <- top_markers_dds(res, n = 5, padj_max = 1e-4, lfc_min = 1)
top5 %>%
    dplyr::mutate(across(-rank, ~ sub(".*\\|", "", .x)))
#> # A tibble: 5 × 9
#>    rank `1`      `2`    `3`    `4`     `5`     `6`    `7`       `8`     
#>   <int> <chr>    <chr>  <chr>  <chr>   <chr>   <chr>  <chr>     <chr>   
#> 1     1 DUT      CCL4L2 KLRG1  SDF2L1  IL7R    LGALS1 LINC00996 ICA1    
#> 2     2 LIG1     CCL4   CST7   HSPA5   NT5E    NA     CD247     NFKBID  
#> 3     3 TMEM106C CCL3L1 DTHD1  TNFRSF9 FLT3LG  NA     MCTP2     MIR155HG
#> 4     4 H2AFX    NA     CD81   TSPAN17 ANXA1   NA     IL2RB     SOD1    
#> 5     5 TMPO     NA     SH2D1A HYOU1   MID1IP1 NA     PTPN12    ICOS
```

The signal is biologically reasonable: cluster 1 is dominated by S-phase
/ DNA-replication genes (`DUT`, `LIG1`, `H2AFX`, `TMPO`), cluster 2 by
effector chemokines (`CCL4`, `CCL3L1`).

## Step 4: pairwise cluster contrasts

`mode = "pairwise"` tests each cluster against each other cluster
individually, returning a long-form table with `group1` and `group2`
columns. This is more conservative than one-vs-all because a true marker
has to differentiate the cluster from *every* other cluster, not just
the average background.

The cost grows quickly — *N* clusters means *N × (N−1)* directional
DESeq2 fits — so it’s worth subsetting to the levels you care about
first. Here we contrast clusters 1, 2, and 3:

``` r

keep <- cc$meta_data$cluster %in% c("1", "2", "3")
md_p  <- cc$meta_data[keep, "cluster", drop = FALSE]
mat_p <- cc$counts_mat[, keep]

res_p <- pseudobulk_deseq2(
    ~cluster, md_p, mat_p,
    verbose = FALSE,
    min_counts_per_sample = 10,
    present_in_min_samples = 5,
    mode = "pairwise"
)
head(res_p)
#>   group1 group2               feature baseMean log2FoldChange lfcSE stat
#> 1      1      3 ENSG00000117632|STMN1    142.2           3.55 0.264 13.5
#> 2      1      3  ENSG00000176890|TYMS     48.3           5.06 0.382 13.2
#> 3      1      2 ENSG00000117632|STMN1    201.2           3.78 0.289 13.1
#> 4      1      2  ENSG00000166508|MCM7     30.7           2.60 0.206 12.6
#> 5      1      2   ENSG00000128951|DUT    142.1           1.65 0.131 12.5
#> 6      1      2 ENSG00000276043|UHRF1     12.7           4.96 0.400 12.4
#>     pvalue     padj
#> 1 3.03e-41 7.42e-38
#> 2 5.09e-40 8.31e-37
#> 3 5.56e-39 5.52e-35
#> 4 2.03e-36 1.01e-32
#> 5 4.28e-36 1.42e-32
#> 6 2.82e-35 7.01e-32
```

[`summarize_dge_pairs()`](https://immunogenomics.github.io/presto/reference/summarize_dge_pairs.md)
collapses the directional pairs to one row per (group, gene). Pass
`"min"` to keep each gene’s worst comparison (the *most conservative*
effect, useful for high-confidence markers) or `"max"` for the best.

``` r

summarize_dge_pairs(res_p, "min") %>%
    head(10) %>%
    dplyr::mutate(feature = sub(".*\\|", "", feature))
#> [1] "min"
#>      group  feature baseMean log2FoldChange lfcSE  stat   pvalue     padj
#>     <char>   <char>    <num>          <num> <num> <num>    <num>    <num>
#>  1:      1    STMN1   201.22           3.78 0.289  13.1 5.56e-39 5.52e-35
#>  2:      1     TYMS    70.03           4.38 0.381  11.5 1.56e-30 1.94e-27
#>  3:      1   TUBA1B   282.72           2.49 0.222  11.3 2.31e-29 2.55e-26
#>  4:      1     PCNA    49.39           2.16 0.198  10.9 6.87e-28 5.68e-25
#>  5:      1     TUBB   309.64           1.75 0.163  10.7 7.74e-27 5.91e-24
#>  6:      1     MCM7    22.69           2.34 0.219  10.7 1.29e-26 6.30e-24
#>  7:      1    UHRF1     9.14           4.73 0.451  10.5 1.11e-25 4.74e-23
#>  8:      1 TMEM106C    34.53           1.88 0.180  10.4 2.55e-25 1.04e-22
#>  9:      1      DUT   102.42           1.49 0.145  10.2 1.27e-24 4.45e-22
#> 10:      1    ASF1B    13.35           3.48 0.343  10.2 3.25e-24 1.06e-21
```

## Step 5: within-cluster Case vs Control

A more interesting biological question than “what defines each cluster”
is “for each cluster, which genes change between irColitis cases and
healthy controls?”. `mode = "within"` answers this. The first variable
in the formula is the **split variable** (the cluster), and the second
is the **contrast variable** (`case`). Within each level of the split
variable, DESeq2 fits the inner formula (here `~case`).

``` r

md_w <- cc$meta_data[, c("cluster", "case")]

res_w <- pseudobulk_deseq2(
    ~cluster + case, md_w, cc$counts_mat,
    verbose = FALSE,
    min_counts_per_sample = 10,
    present_in_min_samples = 5,
    mode = "within"
)
head(res_w)
#>   group                   feature baseMean log2FoldChange lfcSE stat   pvalue
#> 1     1      ENSG00000158517|NCF1    73.62           1.45 0.194 7.50 6.58e-14
#> 2     1      ENSG00000171867|PRNP    30.67           1.80 0.242 7.45 9.11e-14
#> 3     1     ENSG00000162496|DHRS3     8.52           2.62 0.386 6.77 1.29e-11
#> 4     1      ENSG00000185101|ANO9    10.32           1.76 0.261 6.75 1.52e-11
#> 5     1     ENSG00000265972|TXNIP   202.85           1.34 0.199 6.74 1.60e-11
#> 6     1 ENSG00000268804|LINC02132     9.84           2.78 0.414 6.72 1.78e-11
#>       padj
#> 1 3.38e-11
#> 2 4.45e-11
#> 3 4.51e-09
#> 4 4.63e-09
#> 5 4.74e-09
#> 6 4.99e-09
```

Top genes upregulated in cases per cluster:

``` r

top_markers_dds(res_w, n = 5, padj_max = 1e-4, lfc_min = 1) %>%
    dplyr::mutate(across(-rank, ~ sub(".*\\|", "", .x)))
#> # A tibble: 5 × 9
#>    rank `1`   `2`    `3`   `4`   `5`    `6`    `7`     `8`   
#>   <int> <chr> <chr>  <chr> <chr> <chr>  <chr>  <chr>   <chr> 
#> 1     1 NCF1  SPINK2 LAIR1 DHRS7 SPINK2 GNPTAB LDLRAD4 RPS12 
#> 2     2 PRNP  SORBS3 NA    CA10  NA     DHRS7  NA      RPL30 
#> 3     3 DHRS3 KIFC3  NA    NCF1  NA     TXNIP  NA      EEF1B2
#> 4     4 ANO9  GNPTAB NA    KIFC3 NA     FCER1G NA      MT-ND3
#> 5     5 TXNIP CALHM6 NA    MYBL1 NA     TGFBR1 NA      NA
```

Cluster-1 cases up-regulate `NCF1`, `PRNP`, `TXNIP` —
interferon-response and oxidative-stress genes consistent with inflamed
tissue.

## Tips

- **Filter low-count genes early.** `min_counts_per_sample` and
  `present_in_min_samples` are applied before fitting; tighter filters
  speed up the DESeq2 step substantially.
- **Always include batch / channel covariates when possible.** Add them
  to the formula after the contrast variable (e.g. `~cluster + channel`)
  if the design supports it.
- **`collapse_background = TRUE`** (in `mode = "one_vs_all"`) can
  stabilize results when donor representation is uneven across clusters.
- **`mode = "within"` supports both 2-level and ordinal contrasts.**
  Two-level contrasts (factor or character) get a Wald test on the
  level-vs-reference coefficient. Three or more levels are
  integer-encoded and treated as an ordinal trend.

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
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dplyr_1.2.1         Matrix_1.7-4        presto_1.0.0       
#> [4] data.table_1.18.2.1 Rcpp_1.1.1         
#> 
#> loaded via a namespace (and not attached):
#>  [1] SummarizedExperiment_1.40.0 gtable_0.3.6               
#>  [3] xfun_0.56                   bslib_0.10.0               
#>  [5] ggplot2_4.0.2               htmlwidgets_1.6.4          
#>  [7] rhdf5_2.54.1                Biobase_2.70.0             
#>  [9] lattice_0.22-7              rhdf5filters_1.22.0        
#> [11] vctrs_0.7.2                 tools_4.5.2                
#> [13] generics_0.1.4              stats4_4.5.2               
#> [15] parallel_4.5.2              tibble_3.3.1               
#> [17] pkgconfig_2.0.3             RColorBrewer_1.1-3         
#> [19] S7_0.2.1                    desc_1.4.3                 
#> [21] S4Vectors_0.48.0            lifecycle_1.0.5            
#> [23] compiler_4.5.2              farver_2.1.2               
#> [25] textshaping_1.0.4           DESeq2_1.50.2              
#> [27] Seqinfo_1.0.0               codetools_0.2-20           
#> [29] htmltools_0.5.9             sass_0.4.10                
#> [31] yaml_2.3.12                 pillar_1.11.1              
#> [33] pkgdown_2.2.0               jquerylib_0.1.4            
#> [35] tidyr_1.3.2                 BiocParallel_1.44.0        
#> [37] cachem_1.1.0                DelayedArray_0.36.0        
#> [39] abind_1.4-8                 tidyselect_1.2.1           
#> [41] locfit_1.5-9.12             digest_0.6.39              
#> [43] purrr_1.2.1                 fastmap_1.2.0              
#> [45] grid_4.5.2                  cli_3.6.5                  
#> [47] SparseArray_1.10.2          magrittr_2.0.5             
#> [49] S4Arrays_1.10.0             utf8_1.2.6                 
#> [51] dichromat_2.0-0.1           withr_3.0.2                
#> [53] scales_1.4.0                rmarkdown_2.30             
#> [55] XVector_0.50.0              matrixStats_1.5.0          
#> [57] otel_0.2.0                  ragg_1.5.0                 
#> [59] evaluate_1.0.5              knitr_1.51                 
#> [61] GenomicRanges_1.62.1        IRanges_2.44.0             
#> [63] rlang_1.1.7                 glue_1.8.0                 
#> [65] BiocGenerics_0.56.0         jsonlite_2.0.0             
#> [67] R6_2.6.1                    Rhdf5lib_1.32.0            
#> [69] MatrixGenerics_1.22.0       systemfonts_1.3.1          
#> [71] fs_1.6.6
```
