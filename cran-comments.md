## Resubmission

This is a resubmission. The one remaining `\dontrun{}` example,
`load_ircolitis_cd8()`, now carries a comment at the top of the block
explaining why it is not run, as requested: the function downloads a
~121 MB dataset from GEO and requires the Bioconductor package `rhdf5`,
so it is provided for the tutorial on real data and users opt in by
calling it themselves. All other examples are executable (a few in
`\donttest{}`).

Earlier resubmissions in response to the CRAN review also:

* Explained the acronyms in the Description text ("auROC" / "AUC" ->
  "area under the receiver operating characteristic curve").
* Replaced `\dontrun{}` with `\donttest{}` in the `pseudobulk_deseq2()`
  example (guarded by `requireNamespace("DESeq2")`).
* Removed all modification of the global environment (the internal toy
  data generator no longer touches the user's `.Random.seed`).

## Test environments

* local macOS, R 4.5.2
* GitHub Actions: macOS (release), Windows (release),
  Ubuntu (devel, release, oldrel-1)
* win-builder (devel and release)

## R CMD check results

0 ERRORs, 0 WARNINGs, 1 NOTE (New submission).

## Downstream dependencies

There are currently no reverse dependencies on CRAN.

## Notes for reviewers

* The suggested packages DESeq2, rhdf5, SingleCellExperiment,
  SummarizedExperiment, DelayedArray, and BiocStyle are from
  Bioconductor. All uses are guarded with requireNamespace() and are
  optional.
* The vignettes are pre-computed (see vignettes/precompute.R): they use a
  ~121 MB dataset downloaded from GEO, so the shipped .Rmd files are
  static and build offline without network access.
