## Submission

This is the first submission of presto to CRAN.

## Test environments

* local macOS, R 4.5.2
* GitHub Actions: macOS (release), Windows (release),
  Ubuntu (devel, release, oldrel-1)
* win-builder (devel and release)

## R CMD check results

0 ERRORs, 0 WARNINGs.

There is 1 NOTE flagging possibly mis-spelled words in the DESCRIPTION.
These are spelled correctly:

  * auROC        - area under the ROC curve
  * Scalable
  * Wilcoxon
  * Seurat, SingleCellExperiment - software names (single-quoted)

## Downstream dependencies

There are currently no reverse dependencies on CRAN.

## Notes for reviewers

* The suggested packages DESeq2, rhdf5, SingleCellExperiment,
  SummarizedExperiment, and BiocStyle are from Bioconductor. All uses are
  guarded with requireNamespace() and are optional.
* The vignettes are pre-computed (see vignettes/precompute.R): they use a
  ~121 MB dataset downloaded from GEO, so the shipped .Rmd files are static
  and build offline without network access.
