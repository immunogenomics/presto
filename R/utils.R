#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom dplyr %>%
#' @examples
#' x <- 5 %>% sum(10)
#'
#' @usage lhs \%>\% rhs
#' @return return value of rhs function.
NULL


tidy_results <- function(wide_res, features, groups) {
    res <- Reduce(cbind, lapply(wide_res, as.numeric)) %>% data.frame()
    colnames(res) <- names(wide_res)
    res$feature <- rep(features, times = length(groups))
    res$group <- rep(groups, each = length(features))
    ## Plain column names: the .data pronoun is deprecated in tidyselect
    ## contexts such as select() (#20).
    res %>% dplyr::select(
        "feature", "group", "avgExpr", "logFC", "statistic",
        "auc", "pval", "padj", "pct_in", "pct_out"
    )
}


compute_ustat_sparse <- function(grs, group_nnz, group.size, n_obs) {
    ## grs: groups x features rank sums of the stored (shifted) values.
    ## group_nnz: groups x features count of non-zero observations.
    ## Zeros in a feature share the average rank (n_zero + 1) / 2; add their
    ## contribution, then convert the rank sum to the Mann-Whitney U statistic.
    gnz <- group.size - group_nnz
    zero.ranks <- (n_obs - colSums(group_nnz) + 1) / 2
    t(t(gnz) * zero.ranks) + grs - group.size * (group.size + 1) / 2
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

    ## Fully tied features (e.g. all-zero) have zero rank variance, so z is
    ## 0 / 0. Their U statistic is exactly n1n2 / 2 (no separation at all),
    ## so report p = 1 rather than NaN.
    z[!is.finite(z)] <- 0

    pvals <- matrix(2 * pnorm(-abs(as.numeric(z))), ncol = ncol(z))
    return(pvals)
}


#' Shared error for NA values in the data matrix (#25).
#' @noRd
stop_wilcox_na <- function() {
    stop(
        "X contains NA values. Unlike stats::wilcox.test(), which drops ",
        "missing values per observation, wilcoxauc() cannot handle NAs ",
        "and would return silently incorrect results. Remove or impute ",
        "the NA values before calling wilcoxauc().",
        call. = FALSE
    )
}


#' Per-feature Wilcoxon statistics for one in-memory matrix (or one block
#' of a DelayedMatrix). Returns the U statistic, raw group sums, non-zero
#' counts (each ngroups x nfeature), and the per-feature tie groups.
#' @noRd
wilcox_stats_matrix <- function(X, y, grp0, ngroups, group.size, n_obs,
                                nthreads, transposed) {
    if (is(X, "dgCMatrix")) {
        ## A single kernel folds the former five passes (Matrix::t, ranking,
        ## sumGroups on the ranked matrix, and sumGroups / nnzeroGroups on
        ## the original) into one transpose plus one per-feature ranking.
        ## For transposed (observations x features) input, the columns
        ## already are features and no transpose happens at all (#18).
        rr <- cpp_wilcox_stats_dgc(
            X@x, X@p, X@i,
            nfeature = if (transposed) ncol(X) else nrow(X),
            ncell = n_obs,
            grp0, ngroups,
            nthreads = max(1L, as.integer(nthreads)),
            transposed = transposed
        )
        ustat <- compute_ustat_sparse(rr$grs, rr$nnz, group.size, n_obs)
        list(ustat = ustat, group_sums = rr$sums, group_nnz = rr$nnz,
             ties = rr$ties)
    } else {
        if (transposed) {
            ## observations are rows: use the row-wise group reductions and
            ## skip the transpose inside the ranking kernel
            group_sums <- cpp_sumGroups_dense(X, grp0, ngroups)
            group_nnz <- cpp_nnzeroGroups_dense(X, grp0, ngroups)
            rank_res <- cpp_rank_matrix_dense(X, transposed = TRUE)
        } else {
            aux <- cpp_sumGroups_nnz_dense_T(X, grp0, ngroups)
            group_sums <- aux$sums
            group_nnz <- aux$nnz
            rank_res <- rank_matrix(X)
        }
        grs <- sumGroups(rank_res$X_ranked, y)
        ustat <- grs - group.size * (group.size + 1) / 2
        list(ustat = ustat, group_sums = group_sums, group_nnz = group_nnz,
             ties = rank_res$ties)
    }
}


#' Realize one DelayedMatrix block in memory, sparse when the backend
#' advertises sparsity and a sparse coercion is available.
#' @noRd
realize_block <- function(Xb) {
    is_sp <- FALSE
    if (requireNamespace("DelayedArray", quietly = TRUE)) {
        is_sp <- tryCatch(
            DelayedArray::is_sparse(Xb),
            error = function(e) FALSE
        )
    }
    if (is_sp) {
        out <- tryCatch(
            methods::as(Xb, "dgCMatrix"),
            error = function(e) NULL
        )
        if (!is.null(out)) return(out)
    }
    as.matrix(Xb)
}


#' Block-wise Wilcoxon statistics for disk-backed DelayedMatrix input
#' (#26). Splits the features into blocks sized by the
#' `presto.block.elements` option (default 1e7 matrix entries per block,
#' about 80 MB if a block realizes dense), realizes each block in memory,
#' and reuses the in-memory kernels. All statistics are per-feature, so
#' the blocks stitch together exactly.
#' @noRd
wilcox_stats_delayed <- function(X, y, grp0, ngroups, group.size, n_obs,
                                 nthreads, transposed, verbose) {
    nfeat <- if (transposed) ncol(X) else nrow(X)
    block_elems <- getOption("presto.block.elements", 1e7)
    feats_per_block <- max(1L, as.integer(floor(block_elems / n_obs)))
    starts <- seq.int(1L, nfeat, by = feats_per_block)
    if (verbose) {
        message(sprintf(
            "Processing %d features of a %s in %d block(s)",
            nfeat, class(X)[1], length(starts)
        ))
    }
    pieces <- lapply(starts, function(s) {
        idx <- s:min(nfeat, s + feats_per_block - 1L)
        Xb <- if (transposed) {
            X[, idx, drop = FALSE]
        } else {
            X[idx, , drop = FALSE]
        }
        Xb <- realize_block(Xb)
        if (anyNA(Xb)) stop_wilcox_na()
        wilcox_stats_matrix(Xb, y, grp0, ngroups, group.size, n_obs,
                            nthreads, transposed)
    })
    list(
        ustat = do.call(cbind, lapply(pieces, `[[`, "ustat")),
        group_sums = do.call(cbind, lapply(pieces, `[[`, "group_sums")),
        group_nnz = do.call(cbind, lapply(pieces, `[[`, "group_nnz")),
        ties = do.call(c, lapply(pieces, function(p) as.list(p$ties)))
    )
}


#' Column-wise tied ranks of a matrix
#'
#' Ranks the entries of each column independently using the average
#' rank for ties, and returns the per-column tie group sizes needed
#' for the Wilcoxon variance correction. Used internally by
#' [wilcoxauc()] (on a transposed input, so that rows become
#' observations) but exposed as a fast standalone ranking primitive
#' for sparse and dense numeric matrices.
#'
#' @param X Numeric matrix or `dgCMatrix`.
#'
#' @return List with two elements:
#' \itemize{
#'   \item `X_ranked` - matrix with the same shape as `X` containing
#'     per-column tied ranks.
#'   \item `ties` - list of integer vectors, one per column, giving
#'     the sizes of all tie groups encountered in that column. Used
#'     by the Wilcoxon statistic to correct for ties.
#' }
#'
#' @examples
#' set.seed(42)
#' exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
#'                 dimnames = list(paste0("G", 1:25), NULL))
#' rank_res <- rank_matrix(exprs)
#'
#' @seealso [wilcoxauc()]
#'
#' @export
rank_matrix <- function(X) {
    UseMethod("rank_matrix")
}

#' @rdname rank_matrix
#' @export
rank_matrix.dgCMatrix <- function(X) {
    Xr <- X
    ## Force a deep copy of the values slot: the C++ ranker overwrites it
    ## in place, and the caller's matrix must not be modified.
    Xr@x <- X@x + 0
    ties <- cpp_rank_matrix_dgc(Xr@x, Xr@p, nrow(Xr), ncol(Xr))
    return(list(X_ranked = Xr, ties = ties))
}

#' @rdname rank_matrix
#' @export
rank_matrix.matrix <- function(X) {
    cpp_rank_matrix_dense(X)
}

#' Group-wise sum of a matrix along one axis
#'
#' For each unique value of the grouping vector `y`, sums the
#' corresponding rows (or columns) of `X`. Used internally by
#' [wilcoxauc()] and [collapse_counts()], but exposed as a fast
#' group-wise reduction primitive that works on both dense matrices
#' and `dgCMatrix` sparse inputs.
#'
#' @param X Numeric matrix or `dgCMatrix`.
#' @param y Group label vector. Coerced to integer factor codes.
#' @param MARGIN Whether observations are along rows or columns of `X`.
#'   `MARGIN = 2` (default): observations are rows
#'   (`length(y) == nrow(X)`); rows are summed within each group.
#'   `MARGIN = 1`: observations are columns
#'   (`length(y) == ncol(X)`); columns are summed within each group.
#'
#' @return Numeric matrix of shape `n_groups x n_features`. Row order
#'   matches the integer order of `factor(y)`.
#'
#' @examples
#' set.seed(42)
#' exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
#'                 dimnames = list(paste0("G", 1:25), NULL))
#' y <- rep(c("A", "B", "C"), each = 50)
#' sumGroups_res <- sumGroups(exprs, y, 1)
#' sumGroups_res <- sumGroups(t(exprs), y, 2)
#'
#' @seealso [nnzeroGroups()], [wilcoxauc()]
#'
#' @export
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

#' @rdname sumGroups
#' @export
sumGroups.dgCMatrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_sumGroups_dgc_T(X@x, X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1,
                            length(unique(y)))
    } else {
        cpp_sumGroups_dgc(X@x, X@p, X@i, ncol(X), as.integer(y) - 1,
                        length(unique(y)))
    }
}

#' @rdname sumGroups
#' @export
sumGroups.matrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_sumGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))
    } else {
        cpp_sumGroups_dense(X, as.integer(y) - 1, length(unique(y)))
    }
}



#' Group-wise non-zero counts of a matrix along one axis
#'
#' For each unique value of the grouping vector `y`, counts the number
#' of non-zero entries among the corresponding rows (or columns) of
#' `X`. Used internally by [wilcoxauc()] to compute the
#' percent-expressed columns (`pct_in`, `pct_out`), but exposed as a
#' fast group-wise reduction primitive for both dense and `dgCMatrix`
#' inputs.
#'
#' @param X Numeric matrix or `dgCMatrix`.
#' @param y Group label vector. Coerced to integer factor codes.
#' @param MARGIN Whether observations are along rows or columns of `X`.
#'   `MARGIN = 2` (default): observations are rows
#'   (`length(y) == nrow(X)`). `MARGIN = 1`: observations are columns
#'   (`length(y) == ncol(X)`).
#'
#' @return Integer matrix of shape `n_groups x n_features`, where
#'   entry `(g, j)` is the number of observations in group `g` for
#'   which feature `j` is non-zero.
#'
#' @examples
#' set.seed(42)
#' exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
#'                 dimnames = list(paste0("G", 1:25), NULL))
#' y <- rep(c("A", "B", "C"), each = 50)
#' nnz_res <- nnzeroGroups(exprs, y, 1)
#' nnz_res <- nnzeroGroups(t(exprs), y, 2)
#'
#' @seealso [sumGroups()], [wilcoxauc()]
#'
#' @export
nnzeroGroups <- function(X, y, MARGIN = 2) {
    if (MARGIN == 2 & nrow(X) != length(y)) {
        stop("wrong dims")
    } else if (MARGIN == 1 & ncol(X) != length(y)) {
        stop("wrong dims")
    }
    UseMethod("nnzeroGroups")
}

#' @rdname nnzeroGroups
#' @export
nnzeroGroups.dgCMatrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_nnzeroGroups_dgc_T(X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1,
                            length(unique(y)))
    } else {
        cpp_nnzeroGroups_dgc(X@p, X@i, ncol(X), as.integer(y) - 1,
                            length(unique(y)))
    }
}

#' @rdname nnzeroGroups
#' @export
nnzeroGroups.matrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_nnzeroGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))
    } else {
        cpp_nnzeroGroups_dense(X, as.integer(y) - 1, length(unique(y)))
    }
}
