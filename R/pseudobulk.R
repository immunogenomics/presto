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

#' @export
collapse_counts <- function(counts_mat, meta_data, varnames, min_cells_per_group=0) {
    ## give each unique row a hash value for indexing
    hash <- compute_hash(meta_data, varnames)
    idx_keep <- which(!is.na(hash))
    hash <- hash[idx_keep]
    hash <- factor(sprintf('sample_%d', as.integer(hash)))
    meta_data <- meta_data[idx_keep, ]
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
    dplyr::select(-N) %>% 
    data.frame()

    row.names(design_collapsed) <- design_collapsed$sample_id

    ## sum over samples
    counts_collapsed <- presto:::sumGroups(counts_mat, hash, 1) %>% t()
    row.names(counts_collapsed) <- row.names(counts_mat)
    colnames(counts_collapsed) <- levels(hash)

    ## reorder to match design matrix
    counts_collapsed <- counts_collapsed[, design_collapsed$sample_id]
    design_collapsed$sample_id <- NULL
    return(list(counts_mat = counts_collapsed, meta_data = design_collapsed))
}

#' @export
pseudobulk_pairwise <- function(dge_formula, counts_df, meta_data, contrast_var, vals_test, verbose) {
    lapply(vals_test, function(foreground_id) {
        background_ids <- as.character(unique(meta_data$seurat_annotations)) %>% setdiff(foreground_id)
        lapply(background_ids, function(background_id) {
            if (verbose) {
                message(sprintf('%s vs %s', foreground_id, background_id)) 
            }
            suppressMessages({suppressWarnings({
                idx_use <- which(meta_data[[contrast_var]] %in% c(foreground_id, background_id))
                design <- meta_data[idx_use, ]
                
                design[[contrast_var]] <- factor(ifelse(design[[contrast_var]] == foreground_id,
                                                   paste0('cluster_', foreground_id), 
                                                   'background'))
                ## Do DGE with DESeq2
                dds <- DESeqDataSetFromMatrix(
                    countData = counts_df[, idx_use],
                    colData = design,
                    design = dge_formula) %>% 
                    DESeq2::DESeq()

                ## Get results 
                contrast_name <- grep('cluster.*_vs_background', resultsNames(dds), value = TRUE)
                dge_res <- results(dds, name = contrast_name) %>% 
                        data.frame() %>% 
                        tibble::rownames_to_column('feature') %>% 
                        dplyr::arrange(-stat) %>% 
                        dplyr::mutate(group1 = foreground_id, group2 = background_id)

            })})
            return(dge_res)        
        }) %>% 
            purrr::reduce(rbind) 
    }) %>% 
    purrr::reduce(rbind) %>% 
    dplyr::select(group1, group2, feature, dplyr::everything()) %>% 
    dplyr::arrange(group1, -stat)
}

#' @export
pseudobulk_one_vs_all <- function(dge_formula, counts_df, meta_data, contrast_var, vals_test, collapse_background, verbose) {
    Reduce(rbind, lapply(vals_test, function(foreground_id) {
        if (verbose) {
            message(foreground_id)      
        }
        suppressMessages({suppressWarnings({
            ## setup design 
            design <- meta_data            
            design[[contrast_var]] <- factor(ifelse(design[[contrast_var]] == foreground_id,
                                               paste0('cluster_', foreground_id), 
                                               'background'))
            
            ## background clusters should not be treated as independent observations
            if (collapse_background) {
                res <- collapse_counts(counts_df, design, colnames(design))
                design <- res$meta_data
                counts_df <- res$counts_mat                
            }
                        
            ## Do DGE with DESeq2
            dds <- DESeqDataSetFromMatrix(
                countData = counts_df,
                colData = design,
                design = dge_formula) %>% 
                DESeq2::DESeq()

            ## Get results 
            contrast_name <- grep('cluster.*_vs_background', resultsNames(dds), value = TRUE)
            dge_res <- results(dds, name = contrast_name) %>% 
                    data.frame() %>% 
                    tibble::rownames_to_column('feature') %>% 
                    dplyr::arrange(-stat) %>% 
                    dplyr::mutate(group = foreground_id)
        })})
        return(dge_res)
    })) %>% 
    dplyr::select(group, feature, dplyr::everything()) %>% 
    dplyr::arrange(group, -stat)
    
}



#' @export
pseudobulk_deseq2 <- function(dge_formula, meta_data, counts_df, verbose=TRUE, 
                              min_counts_per_sample=10, present_in_min_samples=5, 
                              collapse_background=TRUE, vals_test=NULL,
                              mode=c('one_vs_all', 'pairwise')[1]) {
    message('WARNING: meta_data should only contain pseudobulk identifying variables')
    
    ## filter low expressed genes
    genes_keep <- which(Matrix::rowSums(counts_df >= min_counts_per_sample) >= present_in_min_samples)
    if (verbose) {
        message(sprintf('Filtered out %d genes, analyzing %d genes', nrow(counts_df) - length(genes_keep), length(genes_keep)))
    }
    counts_df <- counts_df[genes_keep, ]
    
    ## assume that the first variable in formula is the main contrast variable
    all_vars <- unlist(strsplit(tail(as.character(dge_formula), 1), split = ' \\+ '))
    if (verbose) {
        message(sprintf('All vars: %s', paste(all_vars, collapse = ', ')))
    }
    contrast_var <- head(all_vars, 1)
    if (verbose) {
        message(sprintf('Contrast var: %s', contrast_var))
    }
    if (is.null(vals_test)) {
        vals_test <- as.character(unique(meta_data[[contrast_var]]))
    } else {
        if (any(!vals_test %in% unique(meta_data[[contrast_var]]))) {
            stop('vals_test must be values in the contrast var')    
        }
    }
    
    if (mode == 'one_vs_all') {
        res <- pseudobulk_one_vs_all(dge_formula, counts_df, meta_data, contrast_var, vals_test, collapse_background, verbose)
        
    } else if (mode == 'pairwise') {
        res <- pseudobulk_pairwise(dge_formula, counts_df, meta_data, contrast_var, vals_test, verbose)
    }
    return (res)
}

#' @export
top_markers_dds <- function(res, n=10, pval_max=1, padj_max=1, lfc_min=1) {
    res %>% 
        dplyr::filter(
            .data$pvalue <= pval_max & 
            .data$padj <= padj_max  &
            log2FoldChange >= lfc_min
        ) %>%
        dplyr::group_by(.data$group) %>%
        dplyr::top_n(n = n, wt = .data$stat) %>% 
        dplyr::mutate(rank = rank(-.data$stat, ties.method = 'random')) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(.data$feature, .data$group, .data$rank) %>% 
        tidyr::spread(.data$group, .data$feature, fill = NA) %>% 
        identity()
}

#' @export
summarize_dge_pairs <- function(dge_res, mode=c('min', 'max')[1]) {
    dge_res <- data.table(dge_res)
    switch(
        mode, 
        min = data.table(dge_res)[, head(.SD[order(stat)], 1), by = .(group1, feature)],
        max = data.table(dge_res)[, head(.SD[order(-stat)], 1), by = .(group1, feature)]
    ) %>% 
        dplyr::select(-group2) %>% 
        dplyr::rename(group = group1) %>% 
        dplyr::arrange(group, -stat)
}