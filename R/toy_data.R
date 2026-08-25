#' Deterministic toy counts shared by [toy_seurat()] and [toy_sce()].
#' 20 genes x 300 cells, two cell types, with a few genes shifted per type
#' so marker tests have signal. Saves and restores the caller's RNG state.
#' @noRd
toy_counts <- function() {
    if (exists(".Random.seed", envir = globalenv())) {
        old_seed <- get(".Random.seed", envir = globalenv())
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
    }
    set.seed(42)
    n_genes <- 20
    n_cells <- 300
    cell_type <- rep(c("jurkat", "t293"), length.out = n_cells)
    counts <- matrix(rpois(n_genes * n_cells, lambda = 2),
                     nrow = n_genes)
    ## give each cell type a handful of upregulated marker genes
    counts[1:3, cell_type == "jurkat"] <-
        counts[1:3, cell_type == "jurkat"] +
        rpois(3 * sum(cell_type == "jurkat"), lambda = 3)
    counts[4:6, cell_type == "t293"] <-
        counts[4:6, cell_type == "t293"] +
        rpois(3 * sum(cell_type == "t293"), lambda = 3)
    dimnames(counts) <- list(
        paste0("G", seq_len(n_genes)),
        paste0("cell", seq_len(n_cells))
    )
    list(
        counts = as(counts, "dgCMatrix"),
        meta_data = data.frame(
            cell_type = cell_type,
            row.names = colnames(counts),
            stringsAsFactors = FALSE
        )
    )
}

#' Toy Seurat object for examples and tests
#'
#' Builds a small deterministic `Seurat` object (20 genes x 300 cells, two
#' cell types in the `cell_type` metadata column) with the installed Seurat
#' version, so no serialized object needs to ship with the package or be
#' updated when Seurat changes its internal class structure. Counts are
#' Poisson-simulated with a few upregulated marker genes per cell type, and
#' log-normalized into the `data` layer. The generator restores the caller's
#' random-number state, so calling it does not perturb reproducibility.
#'
#' Requires the suggested package `Seurat`.
#'
#' @return A `Seurat` object with `counts` and `data` layers and a
#'   `cell_type` metadata column with values `"jurkat"` and `"t293"`.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'     object_seurat <- toy_seurat()
#'     head(wilcoxauc(object_seurat, "cell_type"))
#' }
#'
#' @seealso [toy_sce()], [wilcoxauc()]
#'
#' @export
toy_seurat <- function() {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is required for toy_seurat()")
    }
    d <- toy_counts()
    obj <- suppressWarnings(Seurat::CreateSeuratObject(
        counts = d$counts,
        meta.data = d$meta_data
    ))
    suppressMessages(Seurat::NormalizeData(obj, verbose = FALSE))
}

#' Toy SingleCellExperiment object for examples and tests
#'
#' Builds a small deterministic `SingleCellExperiment` (20 genes x 300
#' cells, two cell types in the `cell_type` colData column) with the
#' installed SingleCellExperiment version, so no serialized object needs to
#' ship with the package. The same simulated counts as [toy_seurat()] are
#' stored in the `counts` assay, with log-normalized values in `logcounts`.
#' The generator restores the caller's random-number state.
#'
#' Requires the suggested package `SingleCellExperiment` (Bioconductor).
#'
#' @return A `SingleCellExperiment` with `counts` and `logcounts` assays
#'   and a `cell_type` colData column with values `"jurkat"` and `"t293"`.
#'
#' @examples
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'     object_sce <- toy_sce()
#'     head(wilcoxauc(object_sce, "cell_type"))
#' }
#'
#' @seealso [toy_seurat()], [wilcoxauc()]
#'
#' @export
toy_sce <- function() {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop(
            "Package 'SingleCellExperiment' is required for toy_sce()"
        )
    }
    d <- toy_counts()
    ## simple library-size log-normalization for the logcounts assay
    sf <- Matrix::colSums(d$counts)
    sf <- sf / mean(sf)
    logcounts <- log1p(Matrix::t(Matrix::t(d$counts) / sf))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = d$counts, logcounts = logcounts),
        colData = d$meta_data
    )
    sce
}
