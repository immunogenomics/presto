# Compute a unique integer hash for each row of a data.frame

Treats the supplied columns as a composite key: rows that agree on every
value in `vars_use` receive the same integer hash, distinct combinations
receive distinct hashes. Used internally by
[`collapse_counts()`](https://immunogenomics.github.io/presto/reference/collapse_counts.md)
to identify pseudobulk groups, but exposed for callers who need the same
indexing in their own code.

## Usage

``` r
compute_hash(data_df, vars_use)
```

## Arguments

- data_df:

  A data.frame.

- vars_use:

  Character vector of column names in `data_df`. Each column is coerced
  to factor before hashing.

## Value

Integer vector of length `nrow(data_df)`.

## See also

[`collapse_counts()`](https://immunogenomics.github.io/presto/reference/collapse_counts.md)
