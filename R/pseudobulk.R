#' @export
model.matrix.full <- function(vars_use, data_df, intercept = 0) {    
    discrete_vars <- c()
    for (varname in vars_use) {
        if ("character" %in% class(data_df[[varname]])) {
            data_df[[varname]] <- factor(data_df[[varname]])
            discrete_vars <- c(discrete_vars, varname)
        } else if ("factor" %in% class(data_df[[varname]])) {
            discrete_vars <- c(discrete_vars, varname)            
        } else {
            stop('no support yet for continuous variables')
        }
    }    

    contrasts_df <- lapply(discrete_vars, function(x) {
        diag(nlevels(data.frame(data_df)[, x, drop = TRUE]))
    })
    names(contrasts_df) <- discrete_vars
    for (varname in discrete_vars) {
        colnames(contrasts_df[[varname]]) <- levels(data_df[[varname]])
        row.names(contrasts_df[[varname]]) <- levels(data_df[[varname]])
        
    }
    res <- model.matrix(as.formula(sprintf("~ %d + %s", intercept, paste(vars_use, collapse = "+"))), 
                 data=data_df, contrasts.arg=contrasts_df)    
        
    return(res)
}

#' @export
collapse_counts <- function(counts_mat, meta_data, varnames) {
    design <- model.matrix.full(varnames, meta_data)

    ## give each unique row a hash value for indexing
    hash <- factor(design %*% matrix(2 ^ seq(0, ncol(design) - 1), ncol = 1))
    hash <- factor(sprintf('sample_%d', as.integer(hash)))
    length(unique(hash))

    ## one hot encoded design matrix, sample level
    design_collapsed <- design %>% 
        as_tibble() %>% 
        cbind(sample_id = hash) %>% 
        unique()
    row.names(design_collapsed) <- design_collapsed$sample_id

    ## in case cells were dropped b/c of NA values
    counts_mat <- counts_mat[, row.names(design)]

    counts_collapsed <- presto:::sumGroups(counts_mat, hash, 1) %>% t()
    row.names(counts_collapsed) <- row.names(counts_mat)
    colnames(counts_collapsed) <- levels(hash)

    ## reorder to match design matrix
    counts_collapsed <- counts_collapsed[, design_collapsed$sample_id]
    
    ## recover the 
    meta_collapsed <- data.frame(row.names = row.names(design_collapsed))
    for (varname in varnames) {
        matching_colnames <- grep(sprintf('^%s', varname), colnames(design_collapsed), value = TRUE)
        vals <- matching_colnames[apply(design_collapsed[, matching_colnames] == 1, 1, which)]
        meta_collapsed[[varname]] <- gsub(varname, '', vals)
    }    
    
    return(list(counts_mat = counts_collapsed, meta_data = meta_collapsed))
}


#' @export
pseudobulk_deseq2 <- function(dge_formula, meta_data, counts_df, verbose=TRUE, 
                   min_counts_per_sample=10, present_in_min_samples=5) {
    ## filter low expressed genes
    genes_keep <- which(rowSums(counts_df >= min_counts_per_sample) >= present_in_min_samples)
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
    Reduce(rbind, lapply(unique(meta_data[[contrast_var]]), function(foreground_id) {
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
            res <- collapse_counts(counts_df, design, all_vars)
            design <- res$meta_data
            counts_df <- res$counts_mat
            
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
    dplyr::select(group, feature, everything())

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