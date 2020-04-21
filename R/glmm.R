#' @export 
regress_out_one_gene <- function(formula_full, formula_reduced, design, y, common_nUMI) {
    tryCatch({
        glmer_res <- lme4::glmer(
            formula_full, cbind(design, y), "poisson", 
            control = lme4::glmerControl(calc.derivs = FALSE, optimizer = "nloptwrap")
        )   

        ## simulate all samples to the same read depth 
        design$logUMI <- log(common_nUMI)
        ypred <- lme4:::predict.merMod(
            object = glmer_res, 
            newdata = design, 
            type = 'response',
            re.form = formula_reduced
        ) 
        return(ypred)        
    }, error = function(e) {
        return(rep())
    })
}


#' @export 
regress_out.presto <- function(obj, formula_full, formula_reduced, features = NULL, do_par = TRUE, common_nUMI=1e6) {
    ## TODO: check that formula_reduced only has RHS 
    
    if (do_par) {
        plan(multicore)
        it_fxn <- furrr::future_map
    }
    else {
        it_fxn <- purrr::map
    }
    if (is.null(features)) {
        features <- rownames(counts_mat)    
    }
    lres <- it_fxn(features, function(feature_use) {
        regress_out_one_gene(
            formula_full = formula_full, 
            formula_reduced = formula_reduced,
            design = obj$meta_data, 
            y = obj$counts_mat[feature_use, ], 
            common_nUMI
        )
    })
    
    obj$corr_mat <- Reduce(Matrix::rbind2, lres)
    rownames(obj$corr_mat) <- features
    colnames(obj$corr_mat) <- colnames(obj$counts_mat)
    return(obj)
}



#' @export 
find_markers_glmm_single_gene <- function(dge_formula, design, y, main_effect, nsim) {
    ## Estimate model 
    tryCatch({
        glmer_res <- lme4::glmer(
            dge_formula, cbind(design, y), 'poisson',
            control = lme4::glmerControl(
                calc.derivs=FALSE, ## 1.1X speedup
                optimizer="nloptwrap" ## 2X speedup
            )
        )    
        
        res <- list()

        ## save the learned variance components
        res$var <- as_tibble(lme4::VarCorr(glmer_res)) %>% 
            dplyr::select(grpvar = grp, vcov) 
        ## save all random and fixed effects
        res$ranef <- lme4::ranef(glmer_res, condVar=TRUE) %>% 
            as.data.frame() %>% 
            dplyr::select(-term)
        res$fixef <- data.frame(beta = lme4::fixef(glmer_res)) %>% 
            tibble::rownames_to_column('effect')

        ## get betas and SDs
        res$dge <- res$ranef %>%  
            subset(grpvar == main_effect) %>% 
            dplyr::select(group = grp, beta=condval) %>%
            dplyr::mutate(group = as.character(group)) ## undo factors 
        res$dge$sigma <- as.numeric(apply(as.data.frame(arm::sim(glmer_res, n.sims=nsim)@ranef[main_effect]), 2, sd))

        return(res)
    }, error = function(e) {
        return(NA)
    })
}


#' @export 
find_markers_glmm <- function(
    dge_formula, counts_mat, meta_data, main_effect, features=NULL, do_par=TRUE, nsim=100
) {
    ## TODO: check that main_effect is in formula
    if (do_par) {
        plan(multicore)
        it_fxn <- furrr::future_map
    } else {
        it_fxn <- purrr::map
    }
    if (is.null(features)) {
        features <- rownames(counts_mat)
    }
    
    ## Get results for each gene
    lres <- it_fxn(features, function(feature_use) {
        find_markers_glmm_single_gene(
            dge_formula = dge_formula, 
            design = meta_data, 
            y = counts_mat[feature_use, ], 
            main_effect, 
            nsim
        )
    })
    names(lres) <- features

    ## in case some glmer models failed
    lres <- lres[!is.na(lres)]    
    
    # Then aggregate the three lists together
    common_el <- purrr::reduce(map(lres, names), intersect)
    res <- map(common_el, function(name) {
        map_dfr(lres, name, .id='feature')
    })
    names(res) <- common_el    
    return(res) 
    
}

