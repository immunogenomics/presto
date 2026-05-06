# Collapse counts based on multiple categorical metadata columns

Collapse counts based on multiple categorical metadata columns

## Usage

``` r
collapse_counts(
  counts_mat,
  meta_data,
  varnames,
  min_cells_per_group = 0,
  keep_n = FALSE,
  how = c("sum", "mean")[1]
)
```

## Arguments

- counts_mat:

  counts matrix where columns represent cells and rows represent
  features

- meta_data:

  data.frame containing cell metadata

- varnames:

  subset of \`meta_data\` column names

- min_cells_per_group:

  minimum cells to keep collapsed group

- keep_n:

  keep or drop the \`N\` column containing the number of cells in each
  group. Default is \`FALSE\`

- how:

  method of collapsing counts from groups. \`sum\` or \`mean\`

## Examples

``` r
m <- matrix(sample.int(8, 100*500, replace=TRUE), nrow=100, ncol=500)
rownames(m) <- paste0("G", 1:100)
colnames(m) <- paste0("C", 1:500)
md1 <- sample(c("a", "b"), 500, replace=TRUE)
md2 <- sample(c("c", "d"), 500, replace=TRUE)
df <- data.frame(md1, md2)
data_collapsed <- collapse_counts(m, df, c("md1", "md2"))
head(data_collapsed$counts_mat)
#>    sample_1 sample_0 sample_2 sample_3
#> G1      627      523      550      573
#> G2      663      498      559      554
#> G3      574      517      575      502
#> G4      605      511      551      520
#> G5      639      548      577      546
#> G6      582      554      596      521
```
