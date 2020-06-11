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
        dplyr::full_join(X1, X2, by = terms_join) %>% 
            dplyr::rowwise() %>% 
            dplyr::mutate(beta = sum(beta.x, beta.y, na.rm = TRUE)) %>% 
            dplyr::select(-beta.x, -beta.y)    
    })
}


.merge_sigmas <- function(X1, X2) {
    terms_join <- setdiff(intersect(colnames(X1), colnames(X2)), 'sigma')
    suppressMessages({
    dplyr::full_join(X1, X2, by = terms_join) %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(sigma = sqrt(sum(sigma.x ^ 2, sigma.y ^ 2, na.rm = TRUE))) %>% 
        dplyr::select(-sigma.x, -sigma.y)
    })
}
