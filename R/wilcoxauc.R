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
                              nthreads = 1, ...) {
    ## Check and possibly correct input values
    if (is(X, "dgeMatrix")) X <- as.matrix(X)
    if (is(X, "data.frame")) X <- as.matrix(X)
    if (is(X, "dgTMatrix")) X <- as(X, "dgCMatrix")
    if (is(X, "TsparseMatrix")) X <- as(X, "dgCMatrix")
    if (ncol(X) != length(y)) stop("number of columns of X does not
                                match length of y")
    if (!is.null(groups_use)) {
        idx_use <- which(y %in% intersect(groups_use, y))
        y <- y[idx_use]
        X <- X[, idx_use]
    }

    y <- factor(y)
    idx_use <- which(!is.na(y))
    if (length(idx_use) < length(y)) {
        y <- y[idx_use]
        X <- X[, idx_use]
        if (verbose)
            message("Removing NA values from labels")
    }

    group.size <- as.numeric(table(y))
    if (length(group.size[group.size > 0]) < 2) {
        stop("Must have at least 2 groups defined.")
    }

    if (is.null(row.names(X))) {
        row.names(X) <- paste0("Feature", seq_len(nrow(X)))
    }

    ## Missing values in X would silently corrupt the ranks (the ranking
    ## code sorts values and cannot drop NAs the way stats::wilcox.test
    ## does), so fail loudly instead of returning wrong numbers. See #25.
    if (anyNA(X)) {
        stop(
            "X contains NA values. Unlike stats::wilcox.test(), which drops ",
            "missing values per observation, wilcoxauc() cannot handle NAs ",
            "and would return silently incorrect results. Remove or impute ",
            "the NA values before calling wilcoxauc()."
        )
    }

    ## Compute primary statistics. Group sizes are computed once here (via
    ## tabulate) and reused throughout instead of re-running table(y) at every
    ## use site.
    ngroups <- nlevels(y)
    group.size <- as.numeric(tabulate(y, ngroups))
    n_obs <- length(y)
    grp0 <- as.integer(y) - 1L
    n1n2 <- group.size * (n_obs - group.size)

    if (is(X, "dgCMatrix")) {
        ## A single kernel folds the former five passes (Matrix::t, ranking,
        ## sumGroups on the ranked matrix, and sumGroups / nnzeroGroups on the
        ## original) into one transpose plus one per-feature ranking, returning
        ## the rank sums, raw sums, non-zero counts, and ties together.
        rr <- cpp_wilcox_stats_dgc(
            X@x, X@p, X@i, nrow(X), ncol(X), grp0, ngroups,
            nthreads = max(1L, as.integer(nthreads))
        )
        group_sums <- rr$sums
        group_nnz <- rr$nnz
        ustat <- compute_ustat_sparse(rr$grs, group_nnz, group.size, n_obs)
        ties <- rr$ties
    } else {
        aux <- cpp_sumGroups_nnz_dense_T(X, grp0, ngroups)
        group_sums <- aux$sums
        group_nnz <- aux$nnz
        rank_res <- rank_matrix(X)
        grs <- sumGroups(rank_res$X_ranked, y)
        ustat <- grs - group.size * (group.size + 1) / 2
        ties <- rank_res$ties
    }

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
    return(tidy_results(res_list, row.names(X), levels(y)))
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
        dplyr::select(.data$feature, .data$group, .data$rank) %>%
        tidyr::spread(.data$group, .data$feature, fill = NA)
}
