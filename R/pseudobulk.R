# Register the bare symbols used inside data.table `[...]` expressions, which
# static analysis (R CMD check) cannot see are column names rather than
# undefined globals. dplyr non-standard evaluation elsewhere in this file uses
# the rlang `.data` pronoun instead, so those names do not need registering.
globalVariables(
    names = c(
        ":=", ".N", ".SD", "sample_id", "N", "stat", "group1", "feature"
    ),
    package = "presto",
    add = TRUE
)

#' Compute a unique integer hash for each row of a data.frame
#'
#' Treats the supplied columns as a composite key: rows that agree on
#' every value in `vars_use` receive the same integer hash, distinct
#' combinations receive distinct hashes. Used internally by
#' [collapse_counts()] to identify pseudobulk groups, but exposed for
#' callers who need the same indexing in their own code.
#'
#' @param data_df A data.frame.
#' @param vars_use Character vector of column names in `data_df`. Each
#'   column is coerced to factor before hashing.
#'
#' @return Integer vector of length `nrow(data_df)`.
#'
#' @seealso [collapse_counts()]
#'
#' @export
compute_hash <- function(data_df, vars_use) {
    base <- 1
    hash <- rep(0, nrow(data_df))
    for (varname in vars_use) {
        vals <- factor(data.frame(data_df)[, varname, drop = TRUE])
        nlevel <- nlevels(vals)
        hash <- hash + (as.integer(vals) - 1) * base
        base <- base * nlevel
    }
    return(hash)
}

#' Collapse a single-cell count matrix into pseudobulks
#'
#' Sums (or averages) the columns of a feature-by-cell count matrix
#' according to one or more cell-metadata columns. The result is a
#' feature-by-pseudobulk matrix where each column pools all cells that
#' share the same combination of metadata values (e.g. all cells from a
#' given donor in a given cluster). The resulting matrix is suitable
#' for bulk-RNA-seq tools such as DESeq2, edgeR, or limma. See the
#' `pseudobulk` vignette for an end-to-end walkthrough.
#'
#' @param counts_mat Counts matrix. Rows are features (genes), columns
#'   are cells. Sparse (`dgCMatrix`) or dense.
#' @param meta_data data.frame of cell metadata. Must have one row per
#'   column of `counts_mat`.
#' @param varnames Character vector of column names in `meta_data` that
#'   together define the pseudobulk grouping. Each unique combination of
#'   values across these columns becomes one pseudobulk sample.
#' @param min_cells_per_group Drop pseudobulks containing fewer than
#'   this many cells. Default `0` (keep all).
#' @param keep_n If `TRUE`, retain the per-pseudobulk cell count as
#'   column `N` in the returned `meta_data`. Default `FALSE`.
#' @param how `"sum"` (default) sums counts across cells in each
#'   pseudobulk; `"mean"` divides by `N` to give the mean expression
#'   per cell. `DESeq2` and other count-based regressions expect sums.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item `counts_mat` - feature-by-pseudobulk numeric matrix.
#'   \item `meta_data` - data.frame with one row per pseudobulk
#'     containing the columns named in `varnames` (and `N` if
#'     `keep_n = TRUE`).
#' }
#'
#' @importFrom data.table data.table
#'
#' @examples
#' m <- matrix(sample.int(8, 100*500, replace=TRUE), nrow=100, ncol=500)
#' rownames(m) <- paste0("G", 1:100)
#' colnames(m) <- paste0("C", 1:500)
#' md1 <- sample(c("a", "b"), 500, replace=TRUE)
#' md2 <- sample(c("c", "d"), 500, replace=TRUE)
#' df <- data.frame(md1, md2)
#' data_collapsed <- collapse_counts(m, df, c("md1", "md2"))
#' head(data_collapsed$counts_mat)
#' head(data_collapsed$meta_data)
#'
#' @seealso [pseudobulk_deseq2()], [compute_hash()]
#'
#' @export
#'
collapse_counts <- function(
    counts_mat,
    meta_data,
    varnames,
    min_cells_per_group = 0,
    keep_n = FALSE,
    how = c("sum", "mean")[1]
) {
    ## give each unique row a hash value for indexing
    hash <- compute_hash(meta_data, varnames)
    idx_keep <- which(!is.na(hash))
    hash <- hash[idx_keep]
    hash <- factor(sprintf("sample_%d", as.integer(hash)))
    ## drop = FALSE so a single-column meta_data stays a data.frame
    meta_data <- data.frame(meta_data)[idx_keep, , drop = FALSE]
    counts_mat <- counts_mat[, idx_keep]

    ## one hot encoded design matrix, sample level
    design_collapsed <- data.frame(meta_data)[, varnames, drop = FALSE] %>%
        cbind(sample_id = hash) %>%
        unique()
    design_collapsed <- data.table(meta_data)[
        , varnames, drop = FALSE, with = FALSE
    ][
        , sample_id := hash
    ][
        , N := .N, by = sample_id
    ][
        N >= min_cells_per_group
    ] %>%
    unique() %>%
    data.frame()

    ## sum over samples
    counts_collapsed <- sumGroups(counts_mat, hash, 1) %>% t()
    row.names(counts_collapsed) <- row.names(counts_mat)
    colnames(counts_collapsed) <- levels(hash)

    ## reorder to match design matrix
    counts_collapsed <- counts_collapsed[, design_collapsed$sample_id]
    design_collapsed$sample_id <- NULL

    if (how == "mean") {
        counts_collapsed <- as.matrix(
            counts_collapsed %*% Matrix::Diagonal(x = 1 / design_collapsed$N)
        )
    }
    if (!keep_n) {
        design_collapsed <- dplyr::select(design_collapsed, -N)
    }
    row.names(design_collapsed) <- design_collapsed$sample_id
    return(list(counts_mat = counts_collapsed, meta_data = design_collapsed))
}

#' Pseudobulk DESeq2: pairwise contrasts
#'
#' For each ordered pair of contrast levels in `vals_test`, fits a
#' DESeq2 model on just those two groups of pseudobulks and returns the
#' foreground-vs-background coefficient. This is more conservative than
#' one-vs-all because a marker has to differentiate the foreground from
#' *every* other group, not just the average background. Cost grows as
#' O(N^2) in the number of levels, so subset to the levels of interest
#' first.
#'
#' Most users should call [pseudobulk_deseq2()] with `mode = "pairwise"`
#' rather than this function directly. The companion
#' [summarize_dge_pairs()] collapses the directional pairs to a single
#' row per (gene, group).
#'
#' @inheritParams pseudobulk_deseq2
#' @param contrast_var Name of the contrast column in `meta_data`
#'   (the first term of `dge_formula`).
#'
#' @return data.frame of DESeq2 results with columns `group1`,
#'   `group2`, `feature`, `baseMean`, `log2FoldChange`, `lfcSE`,
#'   `stat`, `pvalue`, `padj`. Each row reports the test where
#'   pseudobulks of `group1` are the foreground and pseudobulks of
#'   `group2` are the background.
#'
#' @seealso [pseudobulk_deseq2()], [pseudobulk_one_vs_all()],
#'   [pseudobulk_within()], [summarize_dge_pairs()]
#'
#' @importFrom purrr reduce
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr arrange mutate select
#'
#' @export
#'
pseudobulk_pairwise <- function(
    dge_formula,
    counts_df,
    meta_data,
    contrast_var,
    vals_test,
    verbose,
    min_counts_per_sample,
    present_in_min_samples
) {
    requireNamespace("DESeq2")
    lapply(vals_test, function(foreground_id) {
        background_ids <- as.character(
            unique(meta_data[[contrast_var]])
        ) %>% setdiff(foreground_id)
        lapply(background_ids, function(background_id) {
            if (verbose) {
                message(sprintf("%s vs %s", foreground_id, background_id))
            }

            idx_use <- which(
                meta_data[[contrast_var]] %in% c(foreground_id, background_id)
            )
            ## genes could be absent in these two groups
            genes_keep <- which(
                Matrix::rowSums(
                    counts_df[, idx_use] >= min_counts_per_sample
                ) >= present_in_min_samples
            )
            ## filter locally
            counts_df <- counts_df[genes_keep, idx_use]

            suppressMessages({suppressWarnings({
                design <- meta_data[idx_use, , drop = FALSE]
                design[[contrast_var]] <- factor(
                    ifelse(design[[contrast_var]] == foreground_id,
                    paste0("cluster_", foreground_id),
                    "background")
                )
                ## Do DGE with DESeq2
                dds <- DESeq2::DESeqDataSetFromMatrix(
                    countData = counts_df,
                    colData = design,
                    design = dge_formula) %>%
                    DESeq2::DESeq()

                ## Get results
                contrast_name <- grep(
                    "cluster.*_vs_background",
                    DESeq2::resultsNames(dds), value = TRUE
                )
                dge_res <- DESeq2::results(dds, name = contrast_name) %>%
                        data.frame() %>%
                        tibble::rownames_to_column("feature") %>%
                        dplyr::arrange(-.data$stat) %>%
                        dplyr::mutate(
                            group1 = foreground_id,
                            group2 = background_id
                        )

            })})
            return(dge_res)
        }) %>%
            purrr::reduce(rbind)
    }) %>%
    purrr::reduce(rbind) %>%
    dplyr::select("group1", "group2", "feature", dplyr::everything()) %>%
    dplyr::arrange(.data$group1, -.data$stat)
}

#' Pseudobulk DESeq2: one-vs-all contrasts
#'
#' For each level of `contrast_var`, fits a DESeq2 model where that
#' level is the foreground and every other level is pooled as the
#' background. Useful for marker-style differential expression: the
#' return is one set of `(log2FoldChange, padj)` per gene per group,
#' giving a quick view of what's elevated in each group relative to
#' the rest.
#'
#' Most users should call [pseudobulk_deseq2()] with
#' `mode = "one_vs_all"` rather than this function directly; the
#' wrapper handles formula parsing, gene-count filtering, and
#' dispatch. See the `pseudobulk` vignette for a worked example.
#'
#' @inheritParams pseudobulk_deseq2
#' @param contrast_var Name of the contrast column in `meta_data`
#'   (the first term of `dge_formula`).
#'
#' @return data.frame of DESeq2 results with columns `group`,
#'   `feature`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`,
#'   `pvalue`, `padj`. Sorted by `stat` descending within each group.
#'
#' @seealso [pseudobulk_deseq2()], [pseudobulk_pairwise()],
#'   [pseudobulk_within()], [top_markers_dds()]
#'
#' @export
#'
pseudobulk_one_vs_all <- function(
    dge_formula,
    counts_df,
    meta_data,
    contrast_var,
    vals_test,
    collapse_background,
    verbose
) {
    requireNamespace("DESeq2")
    Reduce(rbind, lapply(vals_test, function(foreground_id) {
        if (verbose) {
            message(foreground_id)
        }
        suppressMessages({suppressWarnings({
            ## setup design
            design <- meta_data
            design[[contrast_var]] <- factor(
                ifelse(design[[contrast_var]] == foreground_id,
                paste0("cluster_", foreground_id),
                "background")
            )

            ## background clusters should not be treated as independent obs
            if (collapse_background) {
                res <- collapse_counts(counts_df, design, colnames(design))
                design <- res$meta_data
                counts_df <- res$counts_mat
            }

            ## Do DGE with DESeq2
            dds <- DESeq2::DESeqDataSetFromMatrix(
                countData = counts_df,
                colData = design,
                design = dge_formula) %>%
                DESeq2::DESeq()

            ## Get result
            contrast_name <- grep(
                "cluster.*_vs_background",
                DESeq2::resultsNames(dds),
                value = TRUE
            )
            dge_res <- DESeq2::results(dds, name = contrast_name) %>%
                    data.frame() %>%
                    tibble::rownames_to_column("feature") %>%
                    dplyr::arrange(-.data$stat) %>%
                    dplyr::mutate(group = foreground_id)
        })})
        return(dge_res)
    })) %>%
    dplyr::select("group", "feature", dplyr::everything()) %>%
    dplyr::arrange(.data$group, -.data$stat)

}

#' Pseudobulk DESeq2: within-group contrasts
#'
#' Splits the pseudobulks by `split_var` and, within each split level,
#' fits a DESeq2 model whose contrast variable is the second term of
#' `dge_formula`. Used to test for an effect (e.g. case vs. control)
#' restricted to one cluster at a time, where pooling across clusters
#' would mix biology and batch.
#'
#' Two-level contrasts (factor or character) yield the standard
#' `<var>_<level2>_vs_<level1>` Wald coefficient. Three or more levels
#' are integer-encoded and the fit returns an ordinal trend. Most
#' users should call [pseudobulk_deseq2()] with `mode = "within"`
#' rather than this function directly. See the `pseudobulk` vignette
#' for a worked example using case-vs-control DGE within each cell
#' cluster.
#'
#' @inheritParams pseudobulk_deseq2
#' @param split_var Name of the column in `meta_data` to split on
#'   (the first term of `dge_formula`). Each unique level of
#'   `split_var` produces an independent DESeq2 fit.
#'
#' @return data.frame of DESeq2 results with columns `group` (the
#'   value of `split_var`), `feature`, `baseMean`, `log2FoldChange`,
#'   `lfcSE`, `stat`, `pvalue`, `padj`.
#'
#' @seealso [pseudobulk_deseq2()], [pseudobulk_one_vs_all()],
#'   [pseudobulk_pairwise()]
#'
#' @importFrom stats as.formula
#' @export
#'
pseudobulk_within <- function(
    dge_formula,
    counts_df,
    meta_data,
    split_var,
    vals_test,
    verbose,
    min_counts_per_sample,
    present_in_min_samples
) {
    requireNamespace("DESeq2")
    Reduce(rbind, lapply(vals_test, function(group_test) {
        if (verbose) {
            message(group_test)
        }

        ## remove the split-by-variable (e.g. cluster)
        all_vars <- unlist(
            strsplit(tail(as.character(dge_formula), 1), split = " \\+ ")
        )
        dge_formula <- as.formula(
            paste0("~", paste(tail(all_vars, -1), collapse = "+"))
        )
        contrast_var <- all_vars[[2]]

        suppressMessages({suppressWarnings({
            ## setup design
            idx_use <- which(meta_data[[split_var]] == group_test)
            design <- meta_data[idx_use, , drop = FALSE]

            ## Coerce character to factor so nlevels() reflects the data.
            if (is.character(design[[contrast_var]])) {
                design[[contrast_var]] <- factor(design[[contrast_var]])
            }

            ## Two-level factor: keep as factor for a Wald contrast.
            ## More than two levels: treat as ordinal by integer-encoding.
            if (is.factor(design[[contrast_var]]) &&
                nlevels(design[[contrast_var]]) > 2) {
                design[[contrast_var]] <- as.integer(design[[contrast_var]])
            }

            ## genes could be absent in tested group, filter locally
            genes_keep <- which(
                Matrix::rowSums(
                    counts_df[, idx_use] >= min_counts_per_sample
                ) >= present_in_min_samples
            )
            counts_df <- counts_df[genes_keep, idx_use]

            ## Do DGE with DESeq2
            dds <- DESeq2::DESeqDataSetFromMatrix(
                countData = counts_df,
                colData = design,
                design = dge_formula) %>%
                DESeq2::DESeq()

            ## Build the result coefficient name explicitly. For a 2-level
            ## factor it is "<var>_<level2>_vs_<level1>"; for an
            ## integer-encoded ordinal contrast it is just "<var>".
            final <- design[[contrast_var]]
            if (is.factor(final) && nlevels(final) == 2) {
                contrast_name <- paste0(
                    contrast_var, "_",
                    levels(final)[2], "_vs_", levels(final)[1]
                )
            } else {
                contrast_name <- contrast_var
            }
            dge_res <- DESeq2::results(dds, name = contrast_name) %>%
                    data.frame() %>%
                    tibble::rownames_to_column("feature") %>%
                    dplyr::arrange(-.data$stat) %>%
                    dplyr::mutate(group = group_test)
        })})
        return(dge_res)
    })) %>%
    dplyr::select("group", "feature", dplyr::everything()) %>%
    dplyr::arrange(.data$group, -.data$stat)
}

#' Pseudobulk differential expression with DESeq2
#'
#' Runs `DESeq2::DESeq()` on a feature-by-pseudobulk count matrix to
#' identify genes whose expression differs between groups of
#' pseudobulks. Three test designs are supported, selected by `mode`:
#' \itemize{
#'   \item `"one_vs_all"` (default) - test each level of the contrast
#'     variable against the union of all others. Useful for marker
#'     discovery. Dispatched to [pseudobulk_one_vs_all()].
#'   \item `"pairwise"` - test every ordered pair of levels separately.
#'     Useful for high-confidence markers; combine with
#'     [summarize_dge_pairs()] to keep the most conservative pair.
#'     Dispatched to [pseudobulk_pairwise()].
#'   \item `"within"` - split on the first term of `dge_formula` and
#'     test the second term within each split level. Useful for
#'     condition / case-vs-control comparisons restricted to one
#'     cluster at a time. Dispatched to [pseudobulk_within()].
#' }
#' Pseudobulk inputs are typically produced by [collapse_counts()].
#' Genes with low counts across pseudobulks are filtered out before
#' fitting (controlled by `min_counts_per_sample` and
#' `present_in_min_samples`). See the `pseudobulk` vignette for an
#' end-to-end walkthrough on a real single-cell dataset.
#'
#' @param dge_formula One-sided formula such as `~cluster + donor`.
#'   The first term is treated as the contrast variable in
#'   `"one_vs_all"` and `"pairwise"` modes, or as the split variable
#'   in `"within"` mode (in which case the second term becomes the
#'   contrast variable). Additional terms are kept as covariates in
#'   the DESeq2 design.
#' @param meta_data data.frame of pseudobulk metadata. One row per
#'   pseudobulk; should contain only the variables used in
#'   `dge_formula`.
#' @param counts_df Feature-by-pseudobulk integer count matrix. Rows
#'   are features; columns must align with rows of `meta_data`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @param min_counts_per_sample Minimum count per pseudobulk for a
#'   gene to be considered expressed in that pseudobulk. Default `10`.
#' @param present_in_min_samples Minimum number of pseudobulks in
#'   which a gene must reach `min_counts_per_sample` to be retained.
#'   Default `5`.
#' @param collapse_background Used only when `mode = "one_vs_all"`. If
#'   `TRUE`, background pseudobulks are collapsed across the contrast
#'   variable before each DESeq2 fit, which helps when donor
#'   representation is unbalanced across clusters. Default `TRUE`.
#' @param vals_test Character vector of contrast levels to test. If
#'   `NULL` (default), every level of the contrast variable is tested.
#' @param mode One of `"one_vs_all"` (default), `"pairwise"`, or
#'   `"within"`. See description.
#'
#' @return A long-form data.frame of DESeq2 results. Columns include
#'   the group identifier(s) (`group`, or `group1` / `group2` in
#'   `pairwise` mode), `feature`, and the DESeq2 columns `baseMean`,
#'   `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'     ## 40 genes x 300 cells from 2 clusters across 6 donors
#'     m <- matrix(sample.int(8, 40 * 300, replace = TRUE), nrow = 40)
#'     rownames(m) <- paste0("G", 1:40)
#'     colnames(m) <- paste0("C", 1:300)
#'     meta <- data.frame(
#'         cluster = sample(c("a", "b"), 300, replace = TRUE),
#'         donor = sample(paste0("d", 1:6), 300, replace = TRUE)
#'     )
#'
#'     ## collapse cells into per-(cluster, donor) pseudobulks
#'     dc <- collapse_counts(m, meta, c("cluster", "donor"))
#'
#'     ## test each cluster against the rest
#'     res <- pseudobulk_deseq2(
#'         ~cluster,
#'         dc$meta_data["cluster"],
#'         dc$counts_mat,
#'         verbose = FALSE,
#'         present_in_min_samples = 1,
#'         collapse_background = FALSE,
#'         mode = "one_vs_all"
#'     )
#'     head(res)
#' }
#' }
#'
#' @seealso [collapse_counts()], [top_markers_dds()],
#'   [summarize_dge_pairs()]
#'
#' @export
#'
pseudobulk_deseq2 <- function(
    dge_formula,
    meta_data,
    counts_df,
    verbose = TRUE,
    min_counts_per_sample = 10,
    present_in_min_samples = 5,
    collapse_background = TRUE,
    vals_test = NULL,
    mode = c("one_vs_all", "pairwise", "within")[1]
) {
    requireNamespace("DESeq2")
    if (verbose) {
        message(
            "Note: meta_data should contain only the ",
            "pseudobulk identifying variables used in dge_formula"
        )
    }

    ## filter low expressed genes
    genes_keep <- which(
        Matrix::rowSums(
            counts_df >= min_counts_per_sample
        ) >= present_in_min_samples
    )
    if (verbose) {
        message(
            sprintf(
                "Filtered out %d genes, analyzing %d genes",
                nrow(counts_df) - length(genes_keep),
                length(genes_keep)
            )
        )
    }
    counts_df <- counts_df[genes_keep, ]

    ## assume that the first variable in formula is the main contrast variable
    all_vars <- unlist(
        strsplit(tail(as.character(dge_formula), 1), split = " \\+ ")
    )
    if (verbose) {
        message(sprintf("All vars: %s", paste(all_vars, collapse = ", ")))
    }
    contrast_var <- head(all_vars, 1)
    if (verbose) {
        message(sprintf("Contrast var: %s", contrast_var))
    }
    if (is.null(vals_test)) {
        vals_test <- as.character(unique(meta_data[[contrast_var]]))
    } else {
        if (any(!vals_test %in% unique(meta_data[[contrast_var]]))) {
            stop("vals_test must be values in the contrast var")
        }
    }
    res <- switch(
        mode,
        one_vs_all = pseudobulk_one_vs_all(
            dge_formula, counts_df, meta_data,
            contrast_var, vals_test, collapse_background, verbose),
        pairwise = pseudobulk_pairwise(
            dge_formula, counts_df, meta_data, contrast_var, vals_test, verbose,
            min_counts_per_sample, present_in_min_samples),
        within = pseudobulk_within(
            dge_formula, counts_df, meta_data,
            contrast_var, vals_test, verbose,
            min_counts_per_sample, present_in_min_samples)
    )
    return(res)
}

#' Top markers per group from pseudobulk DESeq2 results
#'
#' Filters and ranks the long-form output of [pseudobulk_deseq2()] to
#' give the most distinguishing features per group. The filtering
#' arguments combine multiplicatively, then the top `n` features per
#' group are kept by descending Wald statistic and pivoted into wide
#' form. Counterpart to [top_markers()] for Wilcoxon-based results.
#'
#' @param res Long-form DESeq2 results from [pseudobulk_deseq2()].
#'   Must have columns `group`, `feature`, `pvalue`, `padj`,
#'   `log2FoldChange`, and `stat`.
#' @param n Number of top features to return per group. Default `10`.
#' @param pval_max Filter features with raw `pvalue > pval_max`.
#'   Default `1` (no filter).
#' @param padj_max Filter features with adjusted `padj > padj_max`.
#'   Default `1` (no filter).
#' @param lfc_min Filter features with `log2FoldChange < lfc_min`.
#'   Default `1` keeps only upregulated features; set to `0` to
#'   include downregulated as well, or `-Inf` to disable.
#'
#' @return tibble in wide form: a `rank` column (1..`n`) and one
#'   column per group containing the gene identifier of the top-ranked
#'   feature at that rank. Cells are `NA` for groups that have fewer
#'   than `n` features passing the filters.
#'
#' @seealso [pseudobulk_deseq2()], [top_markers()]
#'
#' @export
top_markers_dds <- function(
    res,
    n = 10,
    pval_max = 1,
    padj_max = 1,
    lfc_min = 1
) {
    res %>%
        dplyr::filter(
            .data$pvalue <= pval_max &
            .data$padj <= padj_max  &
            .data$log2FoldChange >= lfc_min
        ) %>%
        dplyr::group_by(.data$group) %>%
        dplyr::top_n(n = n, wt = .data$stat) %>%
        dplyr::mutate(rank = rank(-.data$stat, ties.method = "random")) %>%
        dplyr::ungroup() %>%
        dplyr::select("feature", "group", "rank") %>%
        dplyr::arrange(.data$rank) %>%
        tidyr::pivot_wider(
            names_from = "group", values_from = "feature", names_sort = TRUE
        )
}

#' Summarize directional pairwise DGE results
#'
#' Collapses the long-form output of `pseudobulk_deseq2(mode = "pairwise")`
#' to a single row per `(group, gene)` by selecting one comparison per
#' gene-group pair. With `mode = "min"` it keeps the worst comparison
#' (most conservative — the gene must beat *every* other group to
#' have a high statistic). With `mode = "max"` it keeps the best
#' (most permissive — the gene need only beat one).
#'
#' @param dge_res data.frame of pairwise results from
#'   [pseudobulk_deseq2()] with `mode = "pairwise"`. Must have columns
#'   `group1`, `group2`, `feature`, and `stat`.
#' @param mode `"min"` (default) keeps the comparison with the
#'   smallest `stat` for each (group, gene); `"max"` keeps the
#'   largest.
#'
#' @return data.frame with one row per (group, gene), sorted by `stat`
#'   descending within each group. The `group2` column is dropped and
#'   `group1` is renamed to `group`.
#'
#' @seealso [pseudobulk_deseq2()], [pseudobulk_pairwise()]
#'
#' @export
#'
summarize_dge_pairs <- function(dge_res, mode=c("min", "max")[1]) {
    dge_res <- data.table(dge_res)
    switch(
        mode,
        min = data.table(dge_res)[
            , head(.SD[order(stat)], 1), by = list(group1, feature)
        ],
        max = data.table(dge_res)[
            , head(.SD[order(-stat)], 1), by = list(group1, feature)
        ]
    ) %>%
        dplyr::select(-dplyr::any_of("group2")) %>%
        dplyr::rename(group = group1) %>%
        dplyr::arrange(.data$group, -.data$stat)
}
