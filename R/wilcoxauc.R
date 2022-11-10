#' Fast Wilcoxon rank sum test and auROC 
#' 
#' Computes auROC and Wilcoxon p-value based on Gaussian approximation. 
#' Inputs can be 
#' \itemize{
#' \item Dense matrix or data.frame
#' \item Sparse matrix, such as dgCMatrix
#' \item Seurat V3 object
#' \item SingleCellExperiment object
#' }
#' For detailed examples, consult the presto vignette. 
#' 
#' @param X A feature-by-sample matrix, Seurat object, or SingleCellExperiment
#'  object
#' @param y vector of group labels. 
#' @param groups_use (optional) which groups from y vector to test. 
#' @param group_by (Seurat & SCE) name of groups variable ('e.g. Cluster').
#' @param assay (Seurat & SCE) name of feature matrix slot (e.g. 'data' or
#'  'logcounts'). 
#' @param seurat_assay (Seurat) name of Seurat Assay (e.g. 'RNA'). 
#' @param verbose boolean, TRUE for warnings and messages. 
#' @param ... input specific parameters. 
#'
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
    X_matrix <- Seurat::GetAssayData(X, assay = seurat_assay, slot = assay)
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
        logcounts <- SingleCellExperiment::logcounts
        standard_assays <- c(
            "normcounts", "logcounts", "cpm", "tpm",
            "weights", "counts")
        standard_assays <- factor(standard_assays, standard_assays)
        available_assays <- names(SummarizedExperiment::assays(X))
        available_assays <- intersect(standard_assays, available_assays)
        if (length(available_assays) == 0) {
            stop("No assays in SingleCellExperiment object")
        } else {
            assay <- available_assays[1]
        }
    }

    X_matrix <- eval(call(assay, X))
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

tidy_results <- function(wide_res, features, groups) {
    res <- Reduce(cbind, lapply(wide_res, as.numeric)) %>% data.frame()
    colnames(res) <- names(wide_res)
    res$feature <- rep(features, times = length(groups))
    res$group <- rep(groups, each = length(features))
    my_cols <- c("feature", "group",
      "avgExpr", "logFC", "statistic", "auc", "pval", "padj", "pct_in", "pct_out")
    res[,my_cols]
}

compute_ustat <- function(Xr, cols, n1n2, group.size) {
    grs <- sumGroups(Xr, cols)
    if (is(Xr, "dgCMatrix")) {
        gnz <- (group.size - nnzeroGroups(Xr, cols))
        zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
        ustat <- t((t(gnz) * zero.ranks)) + grs -
          group.size * (group.size + 1 ) / 2
    } else {
        ustat <- grs - group.size * (group.size + 1 ) / 2
    }
    return(ustat)
}

compute_pval <- function(ustat, ties, N, n1n2) {
    z <- ustat - .5 * n1n2
    z <- z - sign(z) * .5
    .x1 <- N ^ 3 - N
    .x2 <- 1 / (12 * (N^2 - N))
    rhs <- lapply(ties, function(tvals) {
        (.x1 - sum(tvals ^ 3 - tvals)) * .x2
    }) %>% unlist
    usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
    z <- t(z / usigma)
    pvals <- matrix(2 * pnorm(-abs(as.numeric(z))), ncol = ncol(z))
    return(pvals)
}


# rank_matrix
#
# Utility function to rank columns of matrix
#
# @param X feature by observation matrix.
#
# @examples
#
# data(exprs)
# rank_res <- rank_matrix(exprs)
#
# @return List with 2 items
# \itemize{
# \item X_ranked - matrix of entry ranks
# \item ties - list of tied group sizes
# }
# @noRd
rank_matrix <- function(X) {
    UseMethod("rank_matrix")
}

# @rdname rank_matrix
rank_matrix.dgCMatrix <- function(X) {
    Xr <- Matrix(X, sparse = TRUE)
    ties <- cpp_rank_matrix_dgc(Xr@x, Xr@p, nrow(Xr), ncol(Xr))
    return(list(X_ranked = Xr, ties = ties))
}

# @rdname rank_matrix
rank_matrix.matrix <- function(X) {
    cpp_rank_matrix_dense(X)
}


# sumGroups
#
# Utility function to sum over group labels
#
# @param X matrix
# @param y group labels
# @param MARGIN whether observations are rows (=2) or columns (=1)
#
# @examples
#
# data(exprs)
# data(y)
# sumGroups_res <- sumGroups(exprs, y, 1)
# sumGroups_res <- sumGroups(t(exprs), y, 2)
#
# @return Matrix of groups by features
# @noRd
sumGroups <- function(X, y, MARGIN = 2) {
    if (MARGIN == 2 & nrow(X) != length(y)) {
        stop(
            "nrow(X) != length(y) - the number of rows in the matrix is not
            the same length as group labels"
        )
    } else if (MARGIN == 1 & ncol(X) != length(y)) {
        stop(
            "ncol(X) != length(y) - the number of columns in the matrix is not
             the same length as group labels"
        )
    }
    UseMethod("sumGroups")
}

# @rdname sumGroups
sumGroups.dgCMatrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_sumGroups_dgc_T(X@x, X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1,
                            length(unique(y)))
    } else {
        cpp_sumGroups_dgc(X@x, X@p, X@i, ncol(X), as.integer(y) - 1,
                        length(unique(y)))
    }
}

# @rdname sumGroups
sumGroups.matrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_sumGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))
    } else {
        cpp_sumGroups_dense(X, as.integer(y) - 1, length(unique(y)))
    }
}


# @title nnzeroGroups
#
# Utility function to compute number of zeros-per-feature within group
#
# @param X matrix
# @param y group labels
# @param MARGIN whether observations are rows (=2) or columns (=1)
#
# @examples
#
# data(exprs)
# data(y)
# nnz_res <- nnzeroGroups(exprs, y, 1)
# nnz_res <- nnzeroGroups(t(exprs), y, 2)
#
# @return Matrix of groups by features
# @keywords internal
# @noRd
nnzeroGroups <- function(X, y, MARGIN = 2) {
    if (MARGIN == 2 & nrow(X) != length(y)) {
        stop("wrong dims")
    } else if (MARGIN == 1 & ncol(X) != length(y)) {
        stop("wrong dims")
    }
    UseMethod("nnzeroGroups")
}

# @rdname nnzeroGroups
nnzeroGroups.dgCMatrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_nnzeroGroups_dgc_T(X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1,
                            length(unique(y)))
    } else {
        cpp_nnzeroGroups_dgc(X@p, X@i, ncol(X), as.integer(y) - 1,
                            length(unique(y)))
    }
}

# @rdname nnzeroGroups
nnzeroGroups.matrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_nnzeroGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))
    } else {
        cpp_nnzeroGroups_dense(X, as.integer(y) - 1, length(unique(y)))
    }
}

