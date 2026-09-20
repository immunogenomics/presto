## Resubmission

This is a resubmission. To retire the last `\dontrun{}` entirely, the
`load_ircolitis_cd8()` data downloader has been removed from the package:
it only ever provided a realistic dataset for the vignettes, so it now
lives as a build-time helper (vignettes/ircolitis.R, excluded from the
build) used when the pre-computed vignettes are generated. As a result
there are now no `\dontrun{}` examples, and `rhdf5` / `R.utils` are no
longer needed as suggested packages.

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

* The suggested packages DESeq2, SingleCellExperiment,
  SummarizedExperiment, DelayedArray, and BiocStyle are from
  Bioconductor. All uses are guarded with requireNamespace() and are
  optional.
* The vignettes are pre-computed (see vignettes/precompute.R): they use a
  ~121 MB dataset downloaded from GEO by a build-time helper, so the
  shipped .Rmd files are static and build offline without network access.
