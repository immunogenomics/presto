#' Deterministic toy counts shared by [toy_seurat()] and [toy_sce()].
#' 20 genes x 300 cells, two cell types, with a few genes shifted per type
#' so marker tests have signal.
#'
#' Counts come from a self-contained Park-Miller linear congruential
#' generator rather than R's RNG, so the data are reproducible without
#' calling set.seed() -- which would modify the user's `.Random.seed` in
#' the global environment (disallowed by CRAN policy). Every product stays
#' below 2^53, so the arithmetic is exact and identical on all platforms.
#' @noRd
toy_counts <- function() {
    n_genes <- 20L
    n_cells <- 300L
    cell_type <- rep(c("jurkat", "t293"), length.out = n_cells)

    n <- n_genes * n_cells
    draw <- numeric(n)
    state <- 42
    for (k in seq_len(n)) {
        state <- (16807 * state) %% 2147483647
        draw[k] <- state / 2147483647            # uniform in (0, 1)
    }
    counts <- matrix(as.numeric(floor(draw * 4)), nrow = n_genes)  # 0..3

    ## give each cell type a handful of upregulated marker genes
    jurkat <- cell_type == "jurkat"
    counts[1:3, jurkat] <- counts[1:3, jurkat] + 3
    counts[4:6, !jurkat] <- counts[4:6, !jurkat] + 3

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
