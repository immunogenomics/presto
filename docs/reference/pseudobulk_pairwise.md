# Pseudobulk pairwise

Pseudobulk pairwise

## Usage

``` r
pseudobulk_pairwise(
  dge_formula,
  counts_df,
  meta_data,
  contrast_var,
  vals_test,
  verbose,
  min_counts_per_sample,
  present_in_min_samples
)
```

## Arguments

- dge_formula:

  differential gene expression formula for DESeq2

- counts_df:

  counts matrix

- meta_data:

  data.frame of cell metadata

- contrast_var:

  cell metadata column to use for differential gene expression

- vals_test:

  cell metadata columns

- verbose:

  verbose

- min_counts_per_sample:

  minimum counts per sample to include in differential gene expression

- present_in_min_samples:

  minimum samples with gene counts to include in differential gene
  expression
