# Pseudobulk within

Pseudobulk within

## Usage

``` r
pseudobulk_within(
  dge_formula,
  counts_df,
  meta_data,
  split_var,
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

- split_var:

  \-

- vals_test:

  cell metadata columns

- verbose:

  verbose

- min_counts_per_sample:

  minimum counts per sample to include in differential gene expression

- present_in_min_samples:

  minimum samples with gene counts to include in differential gene
  expression
