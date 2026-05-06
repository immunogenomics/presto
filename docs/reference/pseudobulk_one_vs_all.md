# Pseudobulk one versus all

Pseudobulk one versus all

## Usage

``` r
pseudobulk_one_vs_all(
  dge_formula,
  counts_df,
  meta_data,
  contrast_var,
  vals_test,
  collapse_background,
  verbose
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

- collapse_background:

  collapse background counts according to \`contrast_var\`

- verbose:

  verbose
