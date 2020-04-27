# #' split.presto
# #' 
# #' Split presto object
# #' 
# #' @param obj presto object
# #' @param splitvar split by this variable, defined in presto$meta_data
# #' @param .appendname (bool) whether to append name of obj to daughter object names
# #' @param .sep if .appendname is TRUE, sep character for appending
# #' 
# #' @export 
# split.presto <- function(obj, splitvar, .appendname = TRUE, .sep = '_') {
#     meta_list <- pb$meta_data %>% split(pb$meta_data$Tissue)
#     res <- map(names(meta_list), function(name) {
#         if (.appendname) {
#             newname <- paste(pb$name, name, sep = .sep)
#         } else {
#             newname <- name
#         }

#         list(
#             name = newname,
#             meta_data = meta_list[[name]],
#             counts_mat = pb$counts_mat[, meta_list[[name]]$SampleID]
#         )

#         ## TODO: if glmm has splitby variable, split those results too

#     })

#     names(res) <- names(meta_list)
#     return(res)    
# }

# #' glmm.presto
# #' 
# #' Convenience wrapper to find_markers_glmm
# #' 
# #' @param obj presto object
# #' @param dge_formula lme4 valid formula
# #' @param main_effect variable for which to simulate SD. Results stored in `$dge` (currently, must be random effect)
# #' @param analysis_name
# #' @param features
# #' @param do_par parallel will be done with furrr
# #' @param nsim number simulations to estimate SD of main_effect betas
# #' @param overwrite (bool) overwrite model object
# #' 
# #' @export 
# glmm.presto <- function(obj, dge_formula, main_effect, analysis_name, features=NULL, do_par=TRUE, nsim=100, overwrite=TRUE, ...) {    
#     if (analysis_name %in% names(obj$models) & !overwrite) {
#         stop(sprintf('glmm results \"%s\" already exists. Set overwrite=TRUE to replace', analysis_name))
#     }
#     obj$models[[analysis_name]] <- find_markers_glmm(
#         dge_formula, 
#         obj$counts_mat, 
#         obj$meta_data, 
#         main_effect, 
#         features, 
#         do_par, 
#         nsim
#     )

#     return(obj)
# }



# #' map_dfr_rev
# #' 
# #' @param x
# #' @param y
# #' @param bind_id
# #' 
# #' @export 
# map_dfr_rev <- function(y, x, bind_id) {
#     ## I don't think that purrr explicitly has this function
#     map_dfr(x, y, .id = bind_id)        
# }

# #' combine.model
# #' 
# #' @param model_name
# #' @param models
# #' @param bind_id
# #' 
# #' @export 
# combine.model <- function(model_name, models, bind_id) {
# ## functions to combine lme4 models
#     model_list <- map(models, model_name)
#     model_fields <- purrr::reduce(map(model_list, names), intersect)
#     model_new <- map(model_fields, map_dfr_rev, model_list, bind_id)
#     names(model_new) <- model_fields
#     return(model_new)
# }

# #' combine.models
# #' 
# #' Combine list of DGE models 
# #' 
# #' @param obj_list list of model objects
# #' @param bind_id column name to add to combined data_frames (passed to bind_rows .id parameter)
# #' 
# #' @export 
# combine.models <- function(obj_list, bind_id) {
#     ## get models from objects
#     models <- map(obj_list, 'models')
    
#     ## loop over named models
#     models_common <- purrr::reduce(map(models, names), intersect)
#     res <- map(models_common, combine.model, models, bind_id)
#     names(res) <- models_common
#     return(res)
# }

# #' @export 
# meta_analysis.presto <- function(models, type, key_names) {
#     ## TODO: check that feature and group are the same 
#     ## TODO: instead of feature and group, make generic 
#     ## TODO: check that key_names uniquely define rows
    
#     message('CAUTION: if gene names have underscores, pasrsing will be wrong')
#     res <- switch(
#         type, 
#         'fixedeffects' = { 
#             beta_mat <- bind_rows(models, .id = 'temp_group_var')[, c('beta', 'temp_group_var', key_names)] %>% 
#                 tidyr::spread(temp_group_var, beta, fill=0) %>% 
#                 tidyr::unite(key, all_of(key_names), sep = "_") %>% 
#                 tibble::column_to_rownames('key') %>% as.matrix()


#             sd_mat <- bind_rows(models, .id = 'temp_group_var')[, c('sigma', 'temp_group_var', key_names)] %>% 
#                 tidyr::spread(temp_group_var, sigma, fill=0) %>% 
#                 tidyr::unite(key, all_of(key_names), sep = "_") %>% 
#                 tibble::column_to_rownames('key') %>% as.matrix()


#             tests <- intersect(rownames(beta_mat), rownames(sd_mat))
#             beta_mat <- beta_mat[tests, ]
#             sd_mat <- sd_mat[tests, ]
            
            
# #             tests <- models[[1]] %>% 
# #                 tidyr::unite(key, all_of(key_names), sep = '_') %>% 
# #                 with(key)

# #             beta_mat <- models %>% map('beta') %>% dplyr::bind_cols(testid = tests) %>% 
# #                 tibble::column_to_rownames('testid') %>% 
# #                 as.matrix()

# #             sd_mat <- models %>% map('sigma') %>% dplyr::bind_cols(testid = tests) %>% 
# #                 tibble::column_to_rownames('testid') %>% 
# #                 as.matrix()
            
#             ## summary stats
#             compute_fixed_effects(beta_mat, sd_mat) %>% 
#                 tidyr::separate(test, key_names, sep='_')
            
#         },
#         NULL ## default
#     )
#     if (is.null(res)) {
#         stop(sprintf('meta analysis mode %s is not implemented', type))
#     }
#     return(res)
# }


# #' @export 
# compute_fixed_effects <- function(B, S) {
#     w <- 1 / (S^2)
#     beta_fixed <- rowSums(w * B) / rowSums(w)
#     se_fixed <- sqrt(1 / rowSums(w))

#     Q_stat <- rowSums(w * (B - beta_fixed) ^ 2)
#     I2 <- pmax(0, 100 * ((Q_stat - ncol(B) + 1) / Q_stat), na.rm = FALSE) ## if Q=0, I2 cannot be estimated
#     Q_pval <- pchisq(Q_stat, ncol(B) - 1, lower.tail=FALSE)
#     Q_fdr <- p.adjust(Q_pval, 'BH')
    
    
#     wald_stat <- beta_fixed / se_fixed 
# #     wald_pval <- 2 * (1 - pt(abs(wald_stat), ncol(B) - 1)) ## t-test based p value: what's the right df? 
#     wald_pval <- 2 * (1 - pnorm(abs(wald_stat))) ## z-score based p value
#     wald_fdr <- p.adjust(wald_pval, 'BH')
    
#     res <- tibble(
#         test = rownames(B), 
#         beta = beta_fixed, se = se_fixed, 
#         wald_stat, wald_pval, wald_fdr, 
#         Q_stat, I2, Q_pval, Q_fdr
#     )

#     return(res)
# }



# #' @export 
# combine.presto <- function(
#     obj_list, 
#     obj_name,
#     bind_id,
#     meta_analysis=list(NULL, 'fixedeffects')[[1]], 
#     meta_analysis_model=NULL,
#     meta_analysis_keynames=NULL,
#     combine_models=TRUE,
#     ...
# ) {
#     ## Setup names
#     obj <- list()
#     obj$name <- obj_name
#     obj$meta_analyses <- list()
#     names(obj_list) <- map(obj_list, 'name')
#     if (bind_id %in% colnames(obj_list[[1]]$meta_data)) {
#         obj$meta_data <- map(obj_list, 'meta_data') %>% dplyr::bind_rows(.id=NULL)
#     } else {
#         obj$meta_data <- map(obj_list, 'meta_data') %>% dplyr::bind_rows(.id=bind_id)
#     }
#     obj$counts_mat <- map(obj_list, 'counts_mat') %>% purrr::reduce(Matrix::cbind2)
    
#     ## Bind models
#     ## TODO: implement combining models (with bind_rows)
#     if (combine_models) {
#         obj$models <- combine.models(obj_list, bind_id)
#     }
    
#     ## Do meta analyses
#     if (!is.null(meta_analysis)) {
#         obj$meta_analyses[[meta_analysis_model]] <- obj_list %>% 
#             map('models') %>% 
#             map(meta_analysis_model) %>% 
#             map('dge') %>%
#             meta_analysis.presto(meta_analysis, meta_analysis_keynames)
#     }
#     return(obj)
# }




