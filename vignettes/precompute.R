# Pre-compute the vignettes.
#
# Both vignettes use load_ircolitis_cd8(), which downloads a ~121 MB dataset
# from GEO. CRAN builds must never access the network, so we do NOT evaluate
# that code at build time. Instead we knit the *.Rmd.orig sources here, on a
# machine where the data and Bioconductor packages (rhdf5, DESeq2) are
# available, into static *.Rmd files that ship in the package and build
# offline in seconds.
#
# Re-run this whenever you edit a *.Rmd.orig or change code the vignettes call:
#   R CMD INSTALL .            # install the current source first
#   Rscript vignettes/precompute.R
#
# The generated *.Rmd files are what ship; the *.Rmd.orig sources and this
# script are excluded from the build via .Rbuildignore.

library(knitr)

vign_dir <- "vignettes"
if (!dir.exists(vign_dir) && basename(getwd()) == "vignettes") {
    vign_dir <- "."
}

old <- setwd(vign_dir)
on.exit(setwd(old), add = TRUE)

for (name in c("getting-started", "pseudobulk")) {
    message("Knitting ", name)
    knit(
        input  = paste0(name, ".Rmd.orig"),
        output = paste0(name, ".Rmd")
    )
}

message("Done. Generated static vignettes:")
message("  vignettes/getting-started.Rmd")
message("  vignettes/pseudobulk.Rmd")
