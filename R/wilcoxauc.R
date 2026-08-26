#' Fast Wilcoxon rank-sum test and auROC across groups
#'
#' For every (feature, group) pair, computes the Wilcoxon rank-sum
#' statistic comparing observations in that group against all other
#' observations, and the area under the ROC curve as a measure of
#' separability. P-values come from the standard Gaussian approximation
#' to the U statistic with a tie correction. Returns one row per
#' (feature, group) with effect-size and percent-expressed columns
#' alongside the test statistics.
#'
#' Designed to be fast enough to run on whole-genome × hundred-thousand-
#' cell single-cell matrices in seconds. Sparse `dgCMatrix` inputs are
#' processed without densification. Convenience dispatchers extract the
#' counts matrix and group labels from `Seurat` and
#' `SingleCellExperiment` objects. See the `getting-started` vignette
#' for an end-to-end example on a real dataset.
#'
#' @param X Input data. One of:
#' \itemize{
#'   \item a numeric feature-by-observation matrix or `data.frame`,
#'   \item a sparse `dgCMatrix` of the same shape,
#'   \item a disk-backed `DelayedMatrix` (e.g. HDF5-backed, from the
#'     DelayedArray / HDF5Array packages), which is processed in feature
#'     blocks so the whole matrix never has to be loaded in memory,
#'   \item any other matrix-like class with an `as(., "dgCMatrix")`
#'     coercion method (e.g. BPCells), which is converted up front,
#'   \item a `Seurat` (v3+) object,
#'   \item a `SingleCellExperiment` object.
#' }
#'   `X` must not contain `NA` values. Unlike [stats::wilcox.test()],
#'   which drops missing values per observation, `wilcoxauc()` errors on
#'   `NA` input rather than returning silently incorrect results; remove
#'   or impute missing values first.
#' @param y For matrix input, a character/factor vector of group labels
#'   with length equal to `ncol(X)`. Ignored for `Seurat` /
#'   `SingleCellExperiment` input (use `group_by` instead).
#' @param groups_use Optional character vector restricting the test to
#'   a subset of groups in `y` (or `group_by`). Default `NULL` tests
#'   every group.
#' @param group_by For `Seurat` and `SingleCellExperiment` input, name
#'   of the metadata column that holds the group labels (e.g.
#'   `"cluster"`). For `Seurat`, defaults to `Idents(X)`.
#' @param assay For `Seurat`, the layer name within the selected
#'   assay (e.g. `"data"`, `"counts"`, `"scale.data"`). For
#'   `SingleCellExperiment`, the assay name (e.g. `"logcounts"`,
#'   `"counts"`). Defaults pick a sensible value per input class.
#' @param seurat_assay For `Seurat` input, the name of the assay to
#'   pull from (e.g. `"RNA"`). Default `"RNA"`.
#' @param verbose Logical. Print warnings and informational messages.
#'   Default `TRUE`.
#' @param nthreads Number of threads for the per-feature ranking of sparse
#'   (`dgCMatrix`) input. Default `1` (serial). Values above `1` split the
#'   ranking across threads; the result is identical regardless of the
#'   thread count. Only sparse input is parallelized -- dense matrix and
#'   `data.frame` input are always processed serially. When running under
#'   `R CMD check` or on CRAN, keep this at the default so no more than two
#'   cores are used.
#' @param transposed Set to `TRUE` when your observations (cells,
#'   samples) are in the **rows** of `X` and the features in the
#'   columns -- i.e. `X` is the transpose of the default
#'   features-by-observations layout. The test then runs directly on
#'   that layout without materializing a transposed copy, which saves
#'   time and memory on large matrices. Same convention as the
#'   `transposed` argument of scater's `calculatePCA()` and
#'   `calculateUMAP()`. Only applies to matrix-like input (the `Seurat`
#'   / `SingleCellExperiment` dispatchers always extract
#'   features-by-observations). Default `FALSE`.
#' @param ... Passed to the input-specific method.
#'
#' @examples
#' ## generate a tiny toy dataset
#' set.seed(42)
#' exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
#'                 dimnames = list(paste0("G", 1:25), NULL))
#' y <- rep(c("A", "B", "C"), each = 50)
#'
#' ## on a dense matrix
#' head(wilcoxauc(exprs, y))
#'
#' ## restrict the comparison to a subset of groups
#' head(wilcoxauc(exprs, y, c('A', 'B')))
#'
#' ## on a sparse matrix
#' exprs_sparse <- as(exprs, 'dgCMatrix')
#' head(wilcoxauc(exprs_sparse, y))
#'
#' ## on a Seurat object (>= v3)
#' if (requireNamespace("Seurat", quietly = TRUE) &&
#'     packageVersion("Seurat") >= "3.0") {
#'     object_seurat <- toy_seurat()
#'     head(wilcoxauc(object_seurat, 'cell_type'))
#' }
#'
#' ## on a SingleCellExperiment object
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'     object_sce <- toy_sce()
#'     head(wilcoxauc(object_sce, 'cell_type'))
#' }
#'
#' @return table with the following columns:
#' \itemize{
#' \item \strong{feature} - feature name (e.g. gene name).
#' \item \strong{group} - group name.
#' \item \strong{avgExpr} - mean value of feature in group.
#' \item \strong{logFC} - difference of mean feature values between
#' observations in the group vs out of the group. When the input is
#' log-transformed expression (e.g. Seurat's `"data"` layer or
#' `logcounts`), this difference of means is a log fold change. On raw
#' (untransformed) values it is a plain difference of means, not a fold
#' change.
#' \item \strong{statistic} - Wilcoxon rank sum U statistic.
#' \item \strong{auc} - area under the receiver operator curve.
#' \item \strong{pval} - nominal p value.
#' \item \strong{padj} - Benjamini-Hochberg adjusted p value.
#' \item \strong{pct_in} - Percent of observations in the group with non-zero
#' feature value.
#' \item \strong{pct_out} - Percent of observations out of the group with
#' non-zero feature value.
#' }
#'
#' @seealso [top_markers()] to summarize markers per group;
#'   [pseudobulk_deseq2()] for a count-based pseudobulk alternative.
#'
#' @export
wilcoxauc <- function(X, ...) {
    UseMethod("wilcoxauc")
}

#' @rdname wilcoxauc
#' @export
wilcoxauc.seurat <- function(X, ...) {
    stop("wilcoxauc only implemented for Seurat Version 3, please upgrade to
        run.")
}

#' @rdname wilcoxauc
#' @export
wilcoxauc.Seurat <- function(
    X,
    group_by = NULL,
    assay = "data",
    groups_use = NULL,
    seurat_assay = "RNA",
    ...
) {
    requireNamespace("Seurat")
    X_matrix <- Seurat::GetAssayData(X, assay = seurat_assay, layer = assay)
    if (is.null(group_by)) {
        y <- Seurat::Idents(X)
    } else {
        y <- Seurat::FetchData(X, group_by) %>% unlist %>% as.character()
    }
    wilcoxauc(X_matrix, y, groups_use, ...)
}

#' @rdname wilcoxauc
#' @export
wilcoxauc.SingleCellExperiment <- function(
        X, group_by = NULL, assay = NULL, groups_use = NULL, ...
    ) {
    if (is.null(group_by)) {
        stop("Must specify group_by with SingleCellExperiment")
    } else if (!group_by %in% names(SummarizedExperiment::colData(X))) {
        stop("group_by value is not defined in colData.")
    }
    y <- SummarizedExperiment::colData(X)[[group_by]]

    if (is.null(assay)) {
        standard_assays <- c(
            "normcounts", "logcounts", "cpm", "tpm",
            "weights", "counts")
        available_assays <- intersect(
            standard_assays,
            SummarizedExperiment::assayNames(X)
        )
        if (length(available_assays) == 0) {
            stop("No assays in SingleCellExperiment object")
        } else {
            assay <- available_assays[1]
        }
    }

    X_matrix <- SummarizedExperiment::assay(X, assay)
    wilcoxauc(X_matrix, y, groups_use, ...)
}

#' @rdname wilcoxauc
#' @export
wilcoxauc.default <- function(X, y, groups_use = NULL, verbose = TRUE,
                              nthreads = 1, transposed = FALSE, ...) {
    ## Check and possibly correct input values
    if (is(X, "dgeMatrix")) X <- as.matrix(X)
    if (is(X, "data.frame")) X <- as.matrix(X)
    if (is(X, "dgTMatrix")) X <- as(X, "dgCMatrix")
    if (is(X, "TsparseMatrix")) X <- as(X, "dgCMatrix")
    ## Other matrix-like classes (e.g. BPCells): try a sparse coercion so
    ## that anything with an as(., "dgCMatrix") method just works (#26).
    if (!is.matrix(X) && !is(X, "dgCMatrix") &&
        !inherits(X, "DelayedMatrix")) {
        X_class <- class(X)[1]
        X <- tryCatch(
            as(X, "dgCMatrix"),
            error = function(e) {
                stop(
                    "wilcoxauc() does not know how to handle input of ",
                    "class '", X_class, "'. Convert it to a matrix or ",
                    "dgCMatrix first.",
                    call. = FALSE
                )
            }
        )
    }
    n_obs_dim <- if (transposed) nrow(X) else ncol(X)
    if (n_obs_dim != length(y)) {
        ## If the other dimension matches length(y), the matrix is most
        ## likely in the other orientation: say exactly what to change.
        other_dim <- if (transposed) ncol(X) else nrow(X)
        hint <- ""
        if (other_dim == length(y)) {
            hint <- if (transposed) {
                paste0(
                    "\nX has length(y) columns: if it is the default ",
                    "features x observations layout, drop transposed = TRUE."
                )
            } else {
                paste0(
                    "\nX has length(y) rows: if your matrix is ",
                    "observations x features (samples in rows), call ",
                    "wilcoxauc(X, y, transposed = TRUE)."
                )
            }
        }
        stop(
            "The number of observations in X (",
            if (transposed) "rows" else "columns", " = ", n_obs_dim,
            ") does not match length(y) (", length(y), ").", hint,
            call. = FALSE
        )
    }
    if (!is.null(groups_use)) {
        idx_use <- which(y %in% intersect(groups_use, y))
        y <- y[idx_use]
        X <- if (transposed) X[idx_use, ] else X[, idx_use]
    }

    y <- factor(y)
    idx_use <- which(!is.na(y))
    if (length(idx_use) < length(y)) {
        y <- y[idx_use]
        X <- if (transposed) X[idx_use, ] else X[, idx_use]
        if (verbose)
            message("Removing NA values from labels")
    }

    group.size <- as.numeric(table(y))
    if (length(group.size[group.size > 0]) < 2) {
        stop("Must have at least 2 groups defined.")
    }

    ## Feature names live on rows normally, on columns for transposed input.
    if (transposed) {
        if (is.null(colnames(X))) {
            colnames(X) <- paste0("Feature", seq_len(ncol(X)))
        }
        features <- colnames(X)
    } else {
        if (is.null(rownames(X))) {
            rownames(X) <- paste0("Feature", seq_len(nrow(X)))
        }
        features <- rownames(X)
    }

    ## Missing values in X would silently corrupt the ranks (the ranking
    ## code sorts values and cannot drop NAs the way stats::wilcox.test
    ## does), so fail loudly instead of returning wrong numbers. See #25.
    ## DelayedMatrix input is checked per realized block instead, to avoid
    ## an extra full pass over the on-disk data.
    if (!inherits(X, "DelayedMatrix") && anyNA(X)) {
        stop_wilcox_na()
    }

    ## Compute primary statistics. Group sizes are computed once here (via
    ## tabulate) and reused throughout instead of re-running table(y) at every
    ## use site.
    ngroups <- nlevels(y)
    group.size <- as.numeric(tabulate(y, ngroups))
    n_obs <- length(y)
    grp0 <- as.integer(y) - 1L
    n1n2 <- group.size * (n_obs - group.size)

    if (inherits(X, "DelayedMatrix")) {
        ## Disk-backed input (e.g. HDF5): realize and process feature
        ## blocks so the whole matrix never has to fit in memory (#26).
        st <- wilcox_stats_delayed(
            X, y, grp0, ngroups, group.size, n_obs,
            nthreads, transposed, verbose
        )
    } else {
        st <- wilcox_stats_matrix(
            X, y, grp0, ngroups, group.size, n_obs, nthreads, transposed
        )
    }
    ustat <- st$ustat
    group_sums <- st$group_sums
    group_nnz <- st$group_nnz
    ties <- st$ties

    auc <- t(ustat / n1n2)
    pvals <- compute_pval(ustat, ties, n_obs, n1n2)
    fdr <- apply(pvals, 2, function(x) p.adjust(x, "BH"))

    ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
    group_pct <- sweep(group_nnz, 1, group.size, "/") %>% t()
    group_pct_out <- -group_nnz %>%
        sweep(2, colSums(group_nnz), "+") %>%
        sweep(1, n_obs - group.size, "/") %>% t()
    group_means <- sweep(group_sums, 1, group.size, "/") %>% t()
    cs <- colSums(group_sums)
    lfc <- Reduce(cbind, lapply(seq_len(ngroups), function(g) {
        group_means[, g] - ((cs - group_sums[g, ]) / (n_obs - group.size[g]))
    }))

    res_list <- list(auc = auc,
                pval = pvals,
                padj = fdr,
                pct_in = 100 * group_pct,
                pct_out = 100 * group_pct_out,
                avgExpr = group_means,
                statistic = t(ustat),
                logFC = lfc)
    return(tidy_results(res_list, features, levels(y)))
}


#' Top markers per group from wilcoxauc results
#'
#' Filters and ranks the long-form output of [wilcoxauc()] to give the
#' most distinguishing features per group. The filter arguments combine
#' multiplicatively, then the top `n` features per group are kept by
#' descending `auc` and pivoted into wide form. Counterpart to
#' [top_markers_dds()] for DESeq2-based pseudobulk results.
#'
#' @param res Long-form results table from [wilcoxauc()].
#' @param n Number of top markers to return per group. Default `10`.
#' @param auc_min Drop features with `auc < auc_min`. Default `0`
#'   (no filter); set to `0.5` to keep only features that are positive
#'   markers (more highly expressed in-group than out).
#' @param pval_max Drop features with raw `pval > pval_max`. Default `1`.
#' @param padj_max Drop features with adjusted `padj > padj_max`.
#'   Default `1`.
#' @param pct_in_min Minimum percent (0-100) of in-group observations
#'   with non-zero feature value. Default `0`.
#' @param pct_out_max Maximum percent (0-100) of out-of-group
#'   observations with non-zero feature value. Default `100`.
#'
#' @return tibble in wide form: a `rank` column (1..`n`) and one
#'   column per group containing the feature name of the top-ranked
#'   marker at that rank. Cells are `NA` for groups with fewer than
#'   `n` features that pass the filters.
#'
#' @examples
#' set.seed(42)
#' exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
#'                 dimnames = list(paste0("G", 1:25), NULL))
#' y <- rep(c("A", "B", "C"), each = 50)
#'
#' res <- wilcoxauc(exprs, y)
#'
#' ## top 10 markers per group, restricted to nominally significant,
#' ## up-regulated features (auc > 0.5 means in-group > out-of-group).
#' top_markers(res, n = 10, auc_min = 0.5, pval_max = 0.05)
#'
#' @seealso [wilcoxauc()], [top_markers_dds()]
#'
#' @export
top_markers <- function(res, n = 10, auc_min = 0, pval_max = 1, padj_max = 1,
                        pct_in_min = 0, pct_out_max = 100) {
    res %>%
        dplyr::filter(
            .data$pval <= pval_max &
            .data$padj <= padj_max &
            .data$auc >= auc_min &
            .data$pct_in >= pct_in_min &
            .data$pct_out <= pct_out_max
        ) %>%
        dplyr::group_by(.data$group) %>%
        dplyr::top_n(n = n, wt = .data$auc) %>%
        dplyr::mutate(rank = rank(-.data$auc, ties.method = "random")) %>%
        dplyr::ungroup() %>%
        dplyr::select("feature", "group", "rank") %>%
        dplyr::arrange(.data$rank) %>%
        tidyr::pivot_wider(
            names_from = "group", values_from = "feature", names_sort = TRUE
        )
}
