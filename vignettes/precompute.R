# Pre-compute the vignettes.
#
# Both vignettes load a ~121 MB real dataset from GEO via load_ircolitis_cd8(),
# defined in the build-time helper ircolitis.R (sourced below). That helper is
# NOT part of the package -- it and this script are excluded from the build via
# .Rbuildignore -- so presto ships neither a data downloader nor an rhdf5
# dependency. We knit the *.Rmd.orig sources here, on a machine that has
# network access and the Bioconductor package 'rhdf5', into static *.Rmd files
# that ship in the package and build offline in seconds.
#
# Re-run this whenever you edit a *.Rmd.orig or change code the vignettes call:
#   R CMD INSTALL .            # install the current source first
#   Rscript vignettes/precompute.R
#
# The generated *.Rmd files are what ship.

library(knitr)

vign_dir <- "vignettes"
if (!dir.exists(vign_dir) && basename(getwd()) == "vignettes") {
    vign_dir <- "."
}

old <- setwd(vign_dir)
on.exit(setwd(old), add = TRUE)

# Make load_ircolitis_cd8() available to the vignette chunks. Sourced into the
# environment the chunks are knit in, so the *.Rmd.orig can call it directly.
source("ircolitis.R")

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
