# Pseudobulk DESeq2

Pseudobulk DESeq2

## Usage

``` r
pseudobulk_deseq2(
  dge_formula,
  meta_data,
  counts_df,
  verbose = TRUE,
  min_counts_per_sample = 10,
  present_in_min_samples = 5,
  collapse_background = TRUE,
  vals_test = NULL,
  mode = c("one_vs_all", "pairwise", "within")[1]
)
```

## Arguments

- dge_formula:

  differential gene expression formula for DESeq2

- meta_data:

  data.frame of cell metadata

- counts_df:

  A feature-by-sample matrix

- verbose:

  verbose

- min_counts_per_sample:

  minimum counts per sample to include in differential gene expression

- present_in_min_samples:

  minimum samples with gene counts to include in differential gene
  expression

- collapse_background:

  collapse background. Default is \`TRUE\`

- vals_test:

  cell metadata columns

- mode:

  kind of pseudobulk testing to perform. One of \`one_vs_all\`,
  \`pairwise\`, or \`within\`

## Examples

``` r
if (FALSE) { # \dontrun{
    m <- matrix(sample.int(8, 100*500, replace=TRUE),nrow=100, ncol=500)
    rownames(m) <- paste0("G", 1:100)
    colnames(m) <- paste0("C", 1:500)
    md1 <- sample(c("a", "b"), 500, replace=TRUE)
    md2 <- sample(c("c", "d"), 500, replace=TRUE)
    df <- data.frame(md1, md2)
    data_collapsed <- collapse_counts(m, df, c("md1", "md2"))
    res_mat <- pseudobulk_deseq2(
        ~md1 + md1,
        data_collapsed$meta_data,
        data_collapsed$counts_mat,
        verbose = TRUE,
        present_in_min_samples = 1
    )
    head(res_mat)
} # }
```
