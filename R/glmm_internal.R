.make_tidy_beta <- function(varname, object) {
    res <- cbind(object$betanames_df, object$beta) %>% 
        subset(grpvar == varname) %>% 
        dplyr::select(-grpvar) %>% 
        tidyr::gather(feature, beta, -grp, -term) %>% 
        identity()
    
    if (grepl(':', varname)) {
        ## nested variable
        res <- tidyr::separate(res, grp, unlist(strsplit(varname, ':')), sep = ':')        
    } else {
        ## not nested
        res[[varname]] <- res$grp
        res$grp <- NULL
    } 

    return(res)
}

.make_tidy_sigma <- function(varname, object) {
    res <- cbind(object$betanames_df, object$sigma) %>% 
        subset(grpvar == varname) %>% 
        dplyr::select(-grpvar) %>% 
        tidyr::gather(feature, sigma, -grp, -term) %>% 
        identity()
    
    if (grepl(':', varname)) {
        ## nested variable
        res <- tidyr::separate(res, grp, unlist(strsplit(varname, ':')), sep = ':')        
    } else {
        ## not nested
        res[[varname]] <- res$grp
        res$grp <- NULL
    } 

    return(res)
}


.merge_betas <- function(X1, X2) {
    terms_join <- setdiff(intersect(colnames(X1), colnames(X2)), 'beta')
    suppressMessages({
        res <- dplyr::full_join(X1, X2, by = terms_join) 
        res$beta <- replace_na(res$beta.x, 0) + replace_na(res$beta.y, 0)
        res$beta.x <- NULL
        res$beta.y <- NULL        
    })
    return(res)
}


.merge_sigmas <- function(X1, X2) {
    terms_join <- setdiff(intersect(colnames(X1), colnames(X2)), 'sigma')
    suppressMessages({
        res <- dplyr::full_join(X1, X2, by = terms_join) 
        res$sigma <- sqrt(replace_na(res$sigma.x, 0) ^ 2 + replace_na(res$sigma.y, 0) ^ 2)
        res$sigma.x <- NULL
        res$sigma.y <- NULL        
    })
    return(res)
}
