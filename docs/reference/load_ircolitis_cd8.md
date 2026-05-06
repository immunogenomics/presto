# Download and load the ircolitis tissue CD8 demo dataset

Downloads \`GSE206299_ircolitis-tissue-cd8.h5ad.gz\` from GEO on first
call and caches it under \[tools::R_user_dir()\]. Subsequent calls reuse
the cache. Used to demonstrate presto on a realistic single-cell
dataset.

## Usage

``` r
load_ircolitis_cd8(
  cache_dir = tools::R_user_dir("presto", which = "cache"),
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- cache_dir:

  Directory for the cached \`.h5ad\` file. Defaults to
  \`tools::R_user_dir("presto", which = "cache")\`.

- overwrite:

  If \`TRUE\`, re-download even if the file is cached.

- verbose:

  Print download / parsing progress.

## Value

A list with elements:

- \`counts\` - gene-by-cell \`dgCMatrix\` of raw UMI counts.

- \`obs\` - cell metadata \`data.frame\` (categorical fields decoded to
  character; \`UMAP1\`/\`UMAP2\` and \`PC1..PC24\` kept as numeric).

- \`var\` - gene metadata \`data.frame\`.

## Details

Source: Thomas et al. (Nat. Med. 2024), GSE206299. 25,341 colon-tissue
CD8 T cells from immune checkpoint colitis patients and controls.

Requires the suggested package \`rhdf5\` (Bioconductor).

## Examples

``` r
if (FALSE) { # \dontrun{
  d <- load_ircolitis_cd8()
  dim(d$counts)
  table(d$obs$cluster)
  res <- wilcoxauc(d$counts, d$obs$cluster)
} # }
```
