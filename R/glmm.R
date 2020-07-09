# #' @export 
# nlopt <- function(par, fn, lower, upper, control) {
#     .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
#         opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
#         maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
#     list(par = res$solution,
#          fval = res$objective,
#          conv = if (res$status > 0) 0 else res$status,
#          message = res$message
#     )
# }


#' @export 
glmm_uni <- function(feature, formula, design, response, family, nsim, has_offset) {    
    tryCatch({ 
        ## Fit the model 
        model <- fit_model.presto(formula, design, response[feature, ], family)
        
        ## residuals
        epsilon <- as.numeric(residuals(model, 'response'))
        epsilon_pearson <- as.numeric(residuals(model, 'pearson'))

        ## beta for fixed and meta effects
        beta <- c(
            as.numeric(fixef(model)), ## fixed effects
            as.numeric(as.data.frame(ranef(model))$condval) ## random effects
        )
        
        ## sigma for fixed and random effects
        sigma <- c(
            ## as.matrix below needed b/c lme4 returns dpoMatrix, which doesn't play well 
            sqrt(diag(as.matrix(lme4::vcov.merMod(model)))), ## fixed effects
            as.numeric(apply(as.data.frame(arm::sim(model, n.sims = nsim)@ranef), 2, sd)) ## random effects
        )
        
        if (has_offset) {
            beta <- c(1, beta)
            sigma <- c(0, sigma)
        }

        ## prior distributions for random effects 
        ## CAUTION: Residuals are dropped for lmer without as.data.frame
        prior_sd <- as.numeric(as.data.frame(lme4::VarCorr(model))$sdcor)
        if (isGLMM(model)) {
            ## for GLMMs, use SD of pearson residuals
            prior_sd <- c(prior_sd, sd(residuals(model, 'pearson')))
        }
        
#         rm(model)
#         gc(full = TRUE, verbose = FALSE, reset = TRUE)
        
        return(list(
            status = 0, beta = beta, sigma = sigma, epsilon = epsilon, 
            epsilon_pearson = epsilon_pearson, prior_sd = prior_sd
        ))
    }, error = function(e) {
        return(list(status = 1, error = e))
    })
}

#' @export 
fit_model.presto <- function(formula, design, response, family) {    
    ## initialize model on one feature
    if (family == 'nb') {
        model <- lme4::glmer.nb(
            formula, cbind(design, y = response), 
            control = lme4::glmerControl(calc.derivs = FALSE, optimizer = "nloptwrap")
        )        
    } else if (family == 'poisson') {
        model <- lme4::glmer(
            formula, cbind(design, y = response), 'poisson', 
            control = lme4::glmerControl(calc.derivs = FALSE, optimizer = "nloptwrap")
        )
    } else if (family == 'gaussian') {
        model <- lme4::lmer(
            formula, cbind(design, y = response),
        control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)
        )        
    } else {
        stop(sprintf('(G)LMM family `%s` not supported', family))
    }
    
    return(model)
}


#' @export 
presto.presto <- function(
    formula, 
    design, 
    response, 
    size_varname,
    features=NULL, 
    ncore=1, 
    nsim=100, 
    family='poisson'
) {
    if (is.null(features)) {
        features <- rownames(response)
    }

    ## TODO: make a more rigorous check for this
    if (family %in% c('poisson', 'binomial', 'nb')) {
        message('CAUTION: if using GLMM, make sure your counts are integers!')
    }
    
    ## To make downstream things easier, give exposure variable a dedicated name 
    ## TODO: check that SIZE is valid exposure type variable 
    design$EXPOSURE <- design[[size_varname]]
    fstr <- gsub(size_varname, 'EXPOSURE', as.character(formula))
    formula <- as.formula(sprintf('%s~%s', fstr[[2]], fstr[[3]]))

    ## fit an initial model just to get the names
    model_base <- fit_model.presto(formula, design, response[features[[1]], ], family)
    priornames_df <- as.data.frame(VarCorr(model_base))[, 1:3]
    if (isGLMM(model_base)) {
        ## glmer does not include residuals in VarCorr, lmer does
        priornames_df <- rbind(priornames_df, tibble(grp = 'Residual', var1 = NA, var2 = NA))
    }

    betanames_df <- list(
        tibble(grpvar = names(fixef(model_base)), term = grpvar, grp = grpvar),
        as.data.frame(ranef(model_base), stringsAsFactors = FALSE)[, 1:3]
    ) %>% 
        bind_rows()
    
    ## check if exposure if offset or fixed effect
    has_offset <- !all(map_lgl(model_base@resp$offset, identical, 0))
    if (has_offset) {
        betanames_df <- rbind(
            tibble(grpvar = 'EXPOSURE', term = 'EXPOSURE', grp = 'EXPOSURE'),
            betanames_df
        )
    }

    message('SETUP FUTURES')
    
    ## set up parallel machinery 
    features <- intersect(features, rownames(response))
    if (ncore == 1) {
        future::plan(sequential)
    } else if (ncore %in% c(0, Inf)) {
        ncore <- availableCores()
        future::plan(multiprocess)
    } else {
        ## ncore weirdly not recognized by future
        .ncore <<- ncore
        future::plan(future::multiprocess(workers = .ncore))
        rm(.ncore)
    }
    message('RUN THOSE FUTURES')
    
    lres <- furrr::future_map(features, glmm_uni, formula, design, response, family, nsim, has_offset)
    names(lres) <- features
    lres <- lres[!is.na(lres)]

#     ## chunk into `ncore` futures, one for each core
#     chunk_assign <- sample(rep(seq(ncore), length.out = length(features)))
#     futures_list <- split(features, chunk_assign) %>% map(function(features_chunk) {
#         future(
#             expr = {
#                 lres <- map(
#                     features_chunk, ## input
#                     glmm_uni, ## function
#                     formula, design, response[features_chunk, , drop = FALSE], family, nsim, has_offset ## params
#                 )
#                 names(lres) <- features_chunk
#                 return(lres)
#             }, 
#             lazy = FALSE,
#             globals = FALSE
#         )
#     })    
#     message('COLLECT THOSE FUTURES')
#     lres <- Reduce(append, value(futures_list))
#     lres <- lres[which(map(lres, 'status') == 0)]    
    
    message('AGGREGATE MY OWN RESULTS')
    # Aggregate results 
    common_el <- purrr::reduce(map(lres, names), intersect) %>% setdiff('status')
    res <- map(common_el, function(name) {
        as.matrix(map_dfr(lres, name))
    })
    names(res) <- common_el

    ## remember things names
    res$betanames_df <- betanames_df
    res$priornames_df <- priornames_df
    res$meta_data <- design
    
    if (has_offset) {
        res$design <- list(EXPOSURE = model_base@resp$offset, t(model_base@pp$X), model_base@pp$Zt) %>% 
            purrr::reduce(Matrix::rbind2)
    } else {
        res$design <- list(t(model_base@pp$X), model_base@pp$Zt) %>% 
            purrr::reduce(Matrix::rbind2)
    }
    row.names(res$design) <- res$betanames_df$grp
    
    res$response <- response[names(lres), ]

    ## compute mean genes
    res <- genemeans.presto(res, xpm=1e6)
    
    
    ## flags and stuff
    res$has_offset <- has_offset
    res$family <- family
    res$size_varname <- size_varname
    res$nsim <- nsim
    res$formula <- formula
    
    return(res)
}



#' @export 
correct_counts <- function (object, effects_remove, umi_common, verbose = 0) 
{
    if (missing(umi_common)) {
        umi_common <- mean(log(Matrix::colSums(object$response)))
    }
    effects_remove <- setdiff(effects_remove, "EXPOSURE")
    idx_keep <- object$betanames_df %>% 
        tibble::rowid_to_column("idx") %>% 
#         subset(!grepl(paste(effects_remove, collapse = "|"), grpvar)) %>% 
        subset(!grpvar %in% effects_remove) %>% 
        with(idx)
    if (verbose > 0) {
        message("remove")
        object$betanames_df[-idx_keep, ] %>% with(unique(grpvar)) %>% 
            print()
        message("preserve")
        object$betanames_df[idx_keep, ] %>% with(unique(grpvar)) %>% 
            print()
    }
    
    design_keep <- object$design[idx_keep, ]
    design_keep["EXPOSURE", ] <- umi_common
    effect_keep <- exp(Matrix::crossprod(design_keep, object$beta[idx_keep, 
        ]))
    object$corrected <- Matrix::t(effect_keep + object$epsilon)
    object$corrected <- as(object$corrected, class(object$response)[[1]])
    colnames(object$corrected) <- colnames(object$response)
    row.names(object$corrected) <- row.names(object$response)
    return(object)
}



#' @export 
query.presto <- function(object, feature) {
    return(list(
        effect = cbind(object$betanames_df, beta = object$beta[, feature], sigma = object$sigma[, feature]), 
        prior = cbind(object$priornames_df, sigma = object$prior_sd[, feature]) 
    ))
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
I2.presto <- function(object, effect, effect_levels=NULL, within=NULL, within_levels=NULL) {
    ## TODO: save results, check if results were already computed in cache. Option to force recompute. 
    if (is.null(within)) {
        if (!effect %in% object$priornames_df$grp) {
            stop(sprintf('Model does not contain random effect intercept for term %s', effect))
        } 
        if (is.null(effect_levels)) {
            effect_levels <- subset(object$betanames_df, grpvar == effect)$grp
        }
        idx_use <- object$betanames_df %>% 
            tibble::rowid_to_column('idx') %>% 
            subset(grpvar == effect) %>% 
#             subset(grp %in% effect_levels) %>% 
            with(idx)
        
        I2 <- matrix(compute_I2(t(object$beta[idx_use, ]), t(object$sigma[idx_use, ])), ncol = 1)
        colnames(I2) <- 'I2'
    } else {        
        ## find interaction term effect that matches both 
        ## parse that interaction term into 2 effects
        effect_name <- object$priornames_df$grp[match(
            c(paste(within, effect, sep=':'), paste(effect, within, sep=':')), 
            object$priornames_df$grp
        )]
        
        ## TODO: add functionalty for subsetting on effect and within levels
        
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

    I2 <- as.data.frame(I2)
    rownames(I2) <- colnames(object$beta)
    I2 <- tibble::rownames_to_column(I2, 'feature')
    object$I2 <- I2

    ## TODO: cache result is dictionary 
    
    return(object)
}



#' @export 
toptable.presto <- function(
    object, n=10, max_pval=.05, max_fdr=.05, min_beta=0, max_beta=Inf, rank_by=c('wald', 'beta', 'pval')[1]
) {
    ## TODO: options for tidy vs wide
    ## TODO: options to show stats (tidy only)
    X <- dplyr::inner_join(
        cbind(object$betanames_df, object$beta) %>% 
            subset(grpvar == 'Cluster') %>% 
            dplyr::select(-grpvar, -term) %>%
            tidyr::gather(key, beta, -grp),
        cbind(object$betanames_df, object$sigma) %>% 
            subset(grpvar == 'Cluster') %>% 
            dplyr::select(-grpvar, -term) %>%
            tidyr::gather(key, sigma, -grp),
        by = c('grp', 'key')
    ) %>% 
        dplyr::mutate(
            wald = beta / sigma,
            pval = -2 * pnorm(abs(beta / sigma), log.p = TRUE)
#             pval = 2 * (1 - pnorm(abs(beta / sigma)))
        ) %>% 
        dplyr::mutate(fdr = p.adjust(pval, 'BH')) %>% 
        subset(pval < max_pval & beta > min_beta & beta < max_beta & fdr < max_fdr)

    if (rank_by %in% c('pval')) {
        X$rank_val <- -X[[rank_by]]        
    } else {
        X$rank_val <- X[[rank_by]]
    }
    

    res <- X %>%
        dplyr::group_by(.data$grp) %>% 
        dplyr::top_n(n = n, wt = .data$rank_val) %>% 
        dplyr::mutate(rank = rank(-.data$rank_val, ties.method = "random")) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(.data$key, .data$grp, .data$rank) %>% 
        tidyr::spread(.data$grp, .data$key, fill = NA)

    return(res)    
}


#' @export 
genemeans.presto <- function(object, xpm=1e6) {
    ## Special case of marginalizing variables: 
    ##     take out everything except intercept
    ##     give constant value to exposure 
    b0 <- object$beta[which(object$betanames_df$grpvar == '(Intercept)'), ]
    b1 <- object$beta[which(object$betanames_df$grpvar == 'EXPOSURE'), ]
    log_mu <- b0 + log(xpm)*b1
    
    ## if mean counts < 1, set it to 1
    log_mu <- pmax(log_mu, -log(xpm))
    
    object$log_mu <- data.frame(log_mu) %>% 
        tibble::rownames_to_column('feature') 
    
    return(object)
}


#' @export 
effects.presto <- function(object, effects) {
    effects_available <- unique(object$betanames_df$grpvar)
    if (any(!effects %in% effects_available)) {
        missing_vars <- paste(setdiff(effects, effects_available), collapse = ', ')
        stop(sprintf('missing variables in object: %s', missing_vars))
    }
    
    beta_tidy <- effects %>% purrr::map(.make_tidy_beta, object) %>% purrr::reduce(.merge_betas)
    sigma_tidy <- effects %>% purrr::map(.make_tidy_sigma, object) %>% purrr::reduce(.merge_sigmas)
    suppressMessages({
        object$effect <- dplyr::full_join(beta_tidy, sigma_tidy) %>% 
            dplyr::left_join(object$log_mu) %>% 
        
            dplyr::mutate(
                wald = beta/sigma, 
                pval = -2 * pnorm(abs(beta / sigma), log.p = TRUE)
#                 pval = 2 * (1 - pnorm(abs(beta/sigma)))
            ) %>% 
            dplyr::mutate(fdr = p.adjust(pval, "BH"))
    })    
    return(object)
}


#' @export 
plot_ma.presto <- function(object, max_fdr=.05) {
    object$effect %>% 
        ggplot(aes(log_mu, beta, color = fdr < max_fdr)) + 
            geom_point(alpha = .8, size = .5) + 
            theme_test(base_size = 20) + 
            geom_hline(yintercept = 0, linetype = 2, color = 'black') + 
            labs(x = 'Expected Mean (logCPM)', y = 'Log Fold Change') + 
            scale_color_manual(values = c('grey', 'red')) + 
            guides(color = FALSE) + 
            NULL
    
}

#' @export 
plot_volcano.presto <- function(object, max_fdr=.05) {
    object$effect %>% 
        ggplot(aes(beta, -log10(pval), color = fdr < max_fdr)) + 
            geom_point(alpha = .8, size = .5) + 
            theme_test(base_size = 20) + 
            geom_vline(xintercept = 0, linetype = 2, color = 'black') + 
            labs(y = '-log10 p', x = 'Log Fold Change') + 
            scale_color_manual(values = c('grey', 'red')) + 
            guides(color = FALSE) + 
            NULL
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

