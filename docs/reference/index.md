# Package index

## Wilcoxon and auROC

Single-cell marker discovery with the rank-sum test.

- [`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
  : Fast Wilcoxon rank-sum test and auROC across groups
- [`top_markers()`](https://immunogenomics.github.io/presto/reference/top_markers.md)
  : Top markers per group from wilcoxauc results

## Pseudobulk differential expression

Collapse cells to per-donor pseudobulks and run DESeq2 across them.

- [`collapse_counts()`](https://immunogenomics.github.io/presto/reference/collapse_counts.md)
  : Collapse a single-cell count matrix into pseudobulks
- [`pseudobulk_deseq2()`](https://immunogenomics.github.io/presto/reference/pseudobulk_deseq2.md)
  : Pseudobulk differential expression with DESeq2
- [`top_markers_dds()`](https://immunogenomics.github.io/presto/reference/top_markers_dds.md)
  : Top markers per group from pseudobulk DESeq2 results
- [`summarize_dge_pairs()`](https://immunogenomics.github.io/presto/reference/summarize_dge_pairs.md)
  : Summarize directional pairwise DGE results
- [`pseudobulk_one_vs_all()`](https://immunogenomics.github.io/presto/reference/pseudobulk_one_vs_all.md)
  : Pseudobulk DESeq2: one-vs-all contrasts
- [`pseudobulk_pairwise()`](https://immunogenomics.github.io/presto/reference/pseudobulk_pairwise.md)
  : Pseudobulk DESeq2: pairwise contrasts
- [`pseudobulk_within()`](https://immunogenomics.github.io/presto/reference/pseudobulk_within.md)
  : Pseudobulk DESeq2: within-group contrasts

## Matrix utilities

Group-wise reductions used internally by
[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md).

- [`sumGroups()`](https://immunogenomics.github.io/presto/reference/sumGroups.md)
  : Group-wise sum of a matrix along one axis
- [`nnzeroGroups()`](https://immunogenomics.github.io/presto/reference/nnzeroGroups.md)
  : Group-wise non-zero counts of a matrix along one axis
- [`rank_matrix()`](https://immunogenomics.github.io/presto/reference/rank_matrix.md)
  : Column-wise tied ranks of a matrix
- [`compute_hash()`](https://immunogenomics.github.io/presto/reference/compute_hash.md)
  : Compute a unique integer hash for each row of a data.frame

## Demo data

Example datasets used by the vignettes. The first is a real 25,000-cell
dataset downloaded on demand; the others are tiny shipped Seurat /
SingleCellExperiment objects used by the
[`wilcoxauc()`](https://immunogenomics.github.io/presto/reference/wilcoxauc.md)
dispatch tests and examples.

- [`load_ircolitis_cd8()`](https://immunogenomics.github.io/presto/reference/load_ircolitis_cd8.md)
  : Download and load the ircolitis tissue CD8 demo dataset
- [`object_seurat`](https://immunogenomics.github.io/presto/reference/object_seurat.md)
  : Seurat V3 object with fake data
- [`object_sce`](https://immunogenomics.github.io/presto/reference/object_sce.md)
  : SingleCellExperiment object with fake data
