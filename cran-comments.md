## Resubmission

This is a resubmission. In response to the CRAN review, I have:

* Explained the acronyms in the Description text: "auROC" and "AUC" are
  now spelled out as "area under the receiver operating characteristic
  curve".
* Replaced `\dontrun{}` with `\donttest{}` in the `pseudobulk_deseq2()`
  example (guarded by `requireNamespace("DESeq2")`, as DESeq2 is a
  Bioconductor Suggests).
* Removed all modification of the global environment. The internal toy
  data generator previously called `set.seed()` and saved/restored
  `.Random.seed` in `.GlobalEnv`; it now uses a self-contained
  deterministic generator and never touches the user's RNG state.

One `\dontrun{}` remains, in the `load_ircolitis_cd8()` example. This
function downloads a ~121 MB dataset from GEO over the internet and
requires the Bioconductor Suggests `rhdf5` and `R.utils`, so the example
genuinely cannot be executed during checks (no network access, optional
packages). This matches the documented exception for `\dontrun{}`.

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
