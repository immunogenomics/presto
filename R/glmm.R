#' @export 
nlopt <- function(par, fn, lower, upper, control) {
    .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
        opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
        maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
    list(par = res$solution,
         fval = res$objective,
         conv = if (res$status > 0) 0 else res$status,
         message = res$message
    )
}


#' @export 
glmm_uni <- function(feature, model_base, response, nsim) {    
    tryCatch({        
        model <- lme4::refit(model_base, response[feature, ])
        ## residuals
        epsilon <- as.numeric(residuals(model, 'response'))
        epsilon_pearson <- as.numeric(residuals(model, 'pearson'))

        ## beta for fixed and meta effects
        beta <- c(
            1, ## for offset 
            as.numeric(fixef(model)), ## fixed effects
            as.numeric(as.data.frame(ranef(model))$condval) ## random effects
        )

        ## sigma for fixed and random effects
        sigma <- c(
            0, ## for offset
            sqrt(diag(vcov(model))), ## fixed effects
            as.numeric(apply(as.data.frame(arm::sim(model, n.sims = nsim)@ranef), 2, sd)) ## random effects
        )

        ## prior distributions for random effects 
        ## CAUTION: Residuals are dropped for lmer without as.data.frame
        prior_sd <- as.numeric(as.data.frame(VarCorr(model))$vcov)
#         prior_sd <- as.numeric(VarCorr(model))
        if (isGLMM(model)) {
            ## for GLMMs, use SD of pearson residuals
            prior_sd <- c(prior_sd, sd(residuals(model, 'pearson')))
        }

        return(list(beta = beta, sigma = sigma, epsilon = epsilon, epsilon_pearson = epsilon_pearson, prior_sd = prior_sd))        
    }, error = function(e) {
        print(e)
        return(NA)
    })
}


#' @export 
glmm_multi <- function(formula, design, response, features=NULL, parallel=FALSE, nsim=100, family='poisson') {
    if (parallel) {
        furrr::plan(multicore)
        it_fxn <- furrr::future_map
    } else {
        it_fxn <- purrr::map
    }
    if (is.null(features)) {
        features <- rownames(response)
    }
    
    ## initialize model on one feature
    if (family == 'nb') {
        model_base <- lme4::glmer.nb(
            formula, cbind(design, y = response[features[[1]], ]), 
            control = lme4::glmerControl(calc.derivs = FALSE, optimizer = "nloptwrap")
        )        
    } else if (family == 'poisson') {
        model_base <- lme4::glmer(
            formula, cbind(design, y = response[features[[1]], ]), family, 
            control = lme4::glmerControl(calc.derivs = FALSE, optimizer = "nloptwrap")
        )
    } else if (family == 'gaussian') {
        model_base <- lme4::lmer(
            formula, cbind(design, y = response[features[[1]], ]),
        control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)
        )        
    } else {
        stop(sprintf('(G)LMM family `%s` not supported', family))
    }
    priornames_df <- as.data.frame(VarCorr(model_base))[, 1:3]
    if (isGLMM(model_base)) {
        ## glmer does not include residuals in VarCorr, lmer does
        priornames_df <- rbind(priornames_df, tibble(grp = 'Residual', var1 = NA, var2 = NA))
    }

    betanames_df <- list(
        tibble(grpvar = 'OFFSET', term = 'OFFSET', grp = 'OFFSET'),
        tibble(grpvar = names(fixef(model_base)), term = grpvar, grp = grpvar),
        as.data.frame(ranef(model_base), stringsAsFactors = FALSE)[, 1:3]
    ) %>% 
        bind_rows()

    ## then refit with new responses
    lres <- it_fxn(features, glmm_uni, model_base, response, nsim)
    names(lres) <- features
    lres <- lres[!is.na(lres)] 
    
    # Aggregate results 
    common_el <- purrr::reduce(map(lres, names), intersect)
    res <- map(common_el, function(name) {
        as.matrix(map_dfr(lres, name))
    })
    names(res) <- common_el

    ## remember things names
    res$betanames_df <- betanames_df
    res$priornames_df <- priornames_df
    res$meta_data <- design
    res$design <- list(model_base@resp$offset, t(model_base@pp$X), model_base@pp$Zt) %>% 
        purrr::reduce(Matrix::rbind2)
    res$response <- response[names(lres), ]
    
    return(res)
}



#' @export 
correct_counts <- function(object, effects_remove, umi_common, verbose=0) {
    if (missing(umi_common)) {
        umi_common <- exp(mean(object$meta_data$logUMI))
    }
        
    ## always remove read depth (offset or fixed)
    effects_remove <- union(effects_remove, c('OFFSET', 'logUMI') )
    idx_keep <- object$betanames_df %>% 
        tibble::rowid_to_column('idx') %>% 
        subset(!grepl(paste(effects_remove, collapse = '|'), grpvar)) %>% 
        with(idx)    
    
    if (verbose > 0) {
        message('remove')
        object$betanames_df[-idx_keep, ] %>% with(unique(grpvar)) %>% print()
        message('preserve')
        object$betanames_df[idx_keep, ] %>% with(unique(grpvar)) %>% print()        
    }

    effect_keep <- exp(Matrix::crossprod(object$design[idx_keep, ], object$beta[idx_keep, ]))
    object$corrected <- Matrix::t(umi_common * effect_keep + object$epsilon)
    object$corrected <- as(object$corrected, class(object$response)[[1]])
    
    return(object)
}



#' @export 
compute_I2 <- function(B, S) {
    ## WARNING: only computes I2 across rows right now 
    ## TODO: add option for margin, to compute I2 across rows or colummns
    w <- 1 / (S^2)
    beta_fixed <- rowSums(w * B) / rowSums(w)
    Q_stat <- rowSums(w * (B - beta_fixed) ^ 2)
    I2 <- pmax(0, 100 * ((Q_stat - ncol(B) + 1) / Q_stat), na.rm = FALSE) ## if Q=0, I2 cannot be estimated        
    return(I2)    
}


#' @export 
I2.presto <- function(object, effect, within=NULL) {
    ## TODO: save results, check if results were already computed in cache. Option to force recompute. 
    if (is.null(within)) {
        if (!effect %in% object$priornames_df$grp) {
            stop(sprintf('Model does not contain random effect intercept for term %s', effect))
        }        
        idx_use <- object$betanames_df %>% 
            tibble::rowid_to_column('idx') %>% 
            subset(grpvar == effect) %>% 
            with(idx)

        I2 <- matrix(compute_I2(t(object$beta[idx_use, ]), t(object$sigma[idx_use, ])), ncol = 1)
    } else {        
        ## find interaction term effect that matches both 
        ## parse that interaction term into 2 effects
        effect_name <- object$priornames_df$grp[match(
            c(paste(within, effect, sep=':'), paste(effect, within, sep=':')), 
            object$priornames_df$grp
        )]
        effect_name <- effect_name[!is.na(effect_name)]
        if (length(effect_name) == 0) {
            stop(sprintf('Model does not contain interaction term for %s and %s', effect, within))
        }

        ## Now compute heterogeneity within nested group
        X <- object$betanames_df %>% 
            tibble::rowid_to_column('idx') %>% 
            subset(grpvar == effect_name) %>% 
            tidyr::separate(grp, unlist(strsplit(effect_name, ':')), sep = ':')

        I2 <- split(X$idx, X[[within]]) %>% lapply(function(idx_use) {
            compute_I2(t(object$beta[idx_use, ]), t(object$sigma[idx_use, ]))
        }) %>% 
            bind_rows() %>% 
            as.matrix()
    }

    rownames(I2) <- colnames(object$beta)        
    object$I2 <- I2

    ## TODO: cache result is dictionary 
    
    return(object)
}




# #' @export 
# regress_out_one_gene <- function(formula_full, formula_reduced, design, y, common_nUMI) {
#     tryCatch({
#         glmer_res <- lme4::glmer(
#             formula_full, cbind(design, y), "poisson", 
#             control = lme4::glmerControl(calc.derivs = FALSE, optimizer = "nloptwrap")
#         )   

#         ## simulate all samples to the same read depth 
#         design$logUMI <- log(common_nUMI)
#         ypred <- lme4:::predict.merMod(
#             object = glmer_res, 
#             newdata = design, 
#             type = 'response',
#             re.form = formula_reduced
#         ) 
#         yresid <- residuals(glmer_res, type = 'response') 
#         return(ypred + exp(yresid))
#     }, error = function(e) {
#         return(rep(NA, length(y)))
#     })
# }


# #' @export 
# regress_out.presto <- function(obj, formula_full, formula_reduced, features = NULL, do_par = TRUE, common_nUMI=1e6) {
#     ## TODO: check that formula_reduced only has RHS 
    
#     if (do_par) {
#         plan(multicore)
#         it_fxn <- furrr::future_map
#     }
#     else {
#         it_fxn <- purrr::map
#     }
#     if (is.null(features)) {
#         features <- rownames(obj$counts_mat)    
#     }
#     lres <- it_fxn(features, function(feature_use) {
#         regress_out_one_gene(
#             formula_full = formula_full, 
#             formula_reduced = formula_reduced,
#             design = obj$meta_data, 
#             y = obj$counts_mat[feature_use, ], 
#             common_nUMI
#         )
#     })
    
#     obj$corr_mat <- Reduce(Matrix::rbind2, lres)
#     rownames(obj$corr_mat) <- features
#     colnames(obj$corr_mat) <- colnames(obj$counts_mat)
#     return(obj)
# }



# #' @export 
# find_markers_glmm_single_gene <- function(dge_formula, design, y, main_effect, nsim) {
#     ## Estimate model 
#     tryCatch({
#         glmer_res <- lme4::glmer(
#             dge_formula, cbind(design, y), 'poisson',
#             control = lme4::glmerControl(
#                 calc.derivs=FALSE, ## 1.1X speedup
#                 optimizer="nloptwrap" ## 2X speedup
#             )
#         )    
        
#         res <- list()

#         ## save the learned variance components
#         res$var <- as_tibble(lme4::VarCorr(glmer_res)) %>% 
#             dplyr::select(grpvar = grp, vcov) 
#         ## save all random and fixed effects
#         res$ranef <- lme4::ranef(glmer_res, condVar=TRUE) %>% 
#             as.data.frame() %>% 
#             dplyr::select(-term)
#         res$fixef <- data.frame(beta = lme4::fixef(glmer_res)) %>% 
#             tibble::rownames_to_column('effect')

#         ## get betas and SDs
#         res$dge <- res$ranef %>%  
#             subset(grpvar == main_effect) %>% 
#             dplyr::select(group = grp, beta=condval) %>%
#             dplyr::mutate(group = as.character(group)) ## undo factors 
#         res$dge$sigma <- as.numeric(apply(as.data.frame(arm::sim(glmer_res, n.sims=nsim)@ranef[main_effect]), 2, sd))

#         return(res)
#     }, error = function(e) {
#         return(NA)
#     })
# }


# #' @export 
# find_markers_glmm <- function(
#     dge_formula, counts_mat, meta_data, main_effect, features=NULL, do_par=TRUE, nsim=100
# ) {
#     ## TODO: check that main_effect is in formula
#     if (do_par) {
#         plan(multicore)
#         it_fxn <- furrr::future_map
#     } else {
#         it_fxn <- purrr::map
#     }
#     if (is.null(features)) {
#         features <- rownames(counts_mat)
#     }
    
#     ## Get results for each gene
#     lres <- it_fxn(features, function(feature_use) {
#         find_markers_glmm_single_gene(
#             dge_formula = dge_formula, 
#             design = meta_data, 
#             y = counts_mat[feature_use, ], 
#             main_effect, 
#             nsim
#         )
#     })
#     names(lres) <- features

#     ## in case some glmer models failed
#     lres <- lres[!is.na(lres)] 
    
#     # Then aggregate the three lists together
#     common_el <- purrr::reduce(map(lres, names), intersect)
#     res <- map(common_el, function(name) {
#         map_dfr(lres, name, .id='feature')
#     })
#     names(res) <- common_el    
#     return(res) 
    
# }

