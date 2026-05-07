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
#'     data(object_seurat)
#'     object_seurat <- Seurat::UpdateSeuratObject(object_seurat)
#'     head(wilcoxauc(object_seurat, 'cell_type'))
#' }
#'
#' ## on a SingleCellExperiment object
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'     data(object_sce)
#'     head(wilcoxauc(object_sce, 'cell_type'))
#' }
#'
#' @return table with the following columns:
#' \itemize{
#' \item \strong{feature} - feature name (e.g. gene name).
#' \item \strong{group} - group name.
#' \item \strong{avgExpr} - mean value of feature in group.
#' \item \strong{logFC} - log fold change between observations in group vs out.
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
    wilcoxauc(X_matrix, y, groups_use)
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
    wilcoxauc(X_matrix, y, groups_use)
}

#' @rdname wilcoxauc
#' @export
wilcoxauc.default <- function(X, y, groups_use = NULL, verbose = TRUE, ...) {
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

    ## Compute primary statistics
    group.size <- as.numeric(table(y))
    n1n2 <- group.size * (ncol(X) - group.size)
    if (is(X, "dgCMatrix")) {
        rank_res <- rank_matrix(Matrix::t(X))
    } else {
        rank_res <- rank_matrix(X)
    }

    ustat <- compute_ustat(rank_res$X_ranked, y, n1n2, group.size)
    auc <- t(ustat / n1n2)
    pvals <- compute_pval(ustat, rank_res$ties, ncol(X), n1n2)
    fdr <- apply(pvals, 2, function(x) p.adjust(x, "BH"))

    ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
    group_sums <- sumGroups(X, y, 1)
    group_nnz <- nnzeroGroups(X, y, 1)
    group_pct <- sweep(group_nnz, 1, as.numeric(table(y)), "/") %>% t()
    group_pct_out <- -group_nnz %>%
        sweep(2, colSums(group_nnz) , "+") %>% 
        sweep(1, as.numeric(length(y) - table(y)), "/") %>% t()
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
        group_means[, g] - ((cs - group_sums[g, ]) / (length(y) - gs[g]))
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
