rbind.fill.matrix2 <- function (...) 
{
    matrices <- list(...)
    if (length(matrices) == 0) 
        return()
    if (is.list(matrices[[1]]) && !is.matrix(matrices[[1]])) {
        matrices <- matrices[[1]]
    }
    tmp <- unlist(lapply(matrices, is.factor))
    if (any(tmp)) {
        stop("Input ", paste(which(tmp), collapse = ", "), " is a factor and ", 
            "needs to be converted first to either numeric or character.")
    }
    matrices[] <- lapply(matrices, as.matrix)
    lcols <- lapply(matrices, function(x) plyr::amv_dimnames(x)[[2]])
    cols <- unique(unlist(lcols))
    rows <- unlist(lapply(matrices, nrow))
    nrows <- sum(rows)
    output <- matrix(NA, nrow = nrows, ncol = length(cols))
    colnames(output) <- cols
    pos <- matrix(c(cumsum(rows) - rows + 1, rows), ncol = 2)
    for (i in seq_along(rows)) {
        rng <- seq(pos[i, 1], length.out = pos[i, 2])
        output[rng, lcols[[i]]] <- matrices[[i]]
    }
    rownames(output) <- reduce(map(matrices, rownames), c)
    return(output)
}

make_contrast.presto <- function(object, var_contrast, var_nested=NULL, levels_contrast=NULL, levels_nested=NULL) {
    betanames_df <- object$betanames_df
    if (is.null(levels_contrast)) {
        levels_contrast <- betanames_df %>% 
            subset(grpvar == var_contrast) %>% 
            with(unique(grp))
    }
    
    if (is.null(var_nested)) {
        .x <- betanames_df %>% 
            subset(grpvar == var_contrast & grp %in% levels_contrast)
        terms <- with(.x, as.character(glue::glue('{grpvar}.{grp}.{term}')))
        N <- length(terms)
        res <- (1 + 1 / (N - 1)) * diag(nrow = N) + matrix(1 / (1 - N), nrow = N, ncol = N)
        colnames(res) <- terms
        rownames(res) <- .x$grp
    } else {
        if (is.null(levels_contrast)) {
            levels_nested <- betanames_df %>% 
                subset(grpvar == var_nested) %>% 
                with(unique(grp))
        }
        
        if (glue('{var_nested}:{var_contrast}') %in% betanames_df$grpvar) {
            contrast_design <- betanames_df %>% 
                    subset(grpvar == glue('{var_nested}:{var_contrast}')) %>% 
                    tidyr::separate(grp, c(var_nested, var_contrast), sep = ':', remove = FALSE)            
        } else {
            contrast_design <- betanames_df %>% 
                    subset(grpvar == glue('{var_contrast}:{var_nested}')) %>% 
                    tidyr::separate(grp, c(var_contrast, var_nested), sep = ':', remove = FALSE)            
        }
            
        
         contrast_design <- rbind(
            contrast_design %>%
                dplyr::mutate(covmat_name = as.character(glue::glue('{grpvar}.{grp}.{term}'))) %>% 
                dplyr::select(covmat_name, grp, one_of(var_nested, var_contrast)),
            betanames_df %>% 
                subset(grpvar == var_contrast) %>%
                dplyr::mutate(covmat_name = as.character(glue::glue('{grpvar}.{grp}.{term}'))) %>% 
                tidyr::separate(grp, c(var_contrast, var_nested), sep = ':', remove = FALSE) %>% 
                dplyr::select(covmat_name, grp, one_of(var_nested, var_contrast))

        ) 
        res <- setdiff(unique(contrast_design[[var_nested]]), NA) %>% map(function(level_nested_test) {
            .SD <- contrast_design %>% 
                dplyr::filter(is.na(!!sym(var_nested)) | !!sym(var_nested) == level_nested_test) %>% 
                identity()

            N <- length(unique(.SD[[var_contrast]]))
            C <- t(model_matrix(.SD, as.formula(glue('~0 + {var_contrast}')))) * (1 + 1 / (N - 1))
            C <- C + matrix(1 / (1 - N), nrow = nrow(C), ncol = ncol(C))
            colnames(C) <- .SD$covmat_name
            rownames(C) <- gsub(var_contrast, '', rownames(C))
            rownames(C) <- as.character(glue::glue('{rownames(C)}|{level_nested_test}'))
            return(C)
        }) %>% 
            reduce(rbind.fill.matrix2) %>% 
            replace_na(0)
    }
    return(res)
}
