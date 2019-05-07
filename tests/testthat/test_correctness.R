context('Test correctness of Wilcox results')
library(presto)
data(exprs)
data(y)

test_that('presto::wilcoxauc gives same results as stats::wilcox.test', {
    N <- ncol(exprs)
    D <- nrow(exprs)
    
    ## do Wilcox in stats::wilcow.test
    res_base_r <- split(seq_len(N), y) %>% lapply(function(idx) {
        res <- Reduce(rbind, apply(exprs, 1, function(g) {
            wilcox.test(g[idx], g[setdiff(seq_len(N), idx)], exact = FALSE, 
                        correct = TRUE,
                        alternative = 'two.sided') %>% 
                broom::tidy()
        })) %>% data.frame()
        row.names(res) <- row.names(exprs)
        return(res)
    })
    res_base_r <- Reduce(rbind, lapply(names(res_base_r), function(group) {
        res_base_r[[group]] %>% data.frame() %>% 
            tibble::rownames_to_column('feature') %>% 
            dplyr::mutate(group = group)
    })) %>% 
        dplyr::select(feature, group, statistic, p.value)
    
    ## do Wilcox in presto
    res_presto <- wilcoxauc(exprs, y) %>% 
        dplyr::select(feature, group, statistic, p.value = pval)
    
    ## compare U-statistic and p values
    res_joint <- dplyr::inner_join(
        res_base_r, res_presto, 
        by = c('feature', 'group'), suffix = c('_base', '_presto')
    )
    expect_lt(max(abs(res_joint$p.value_base - res_joint$p.value_presto)), 1e-3)
    expect_lt(max(abs(res_joint$statistic_base - res_joint$statistic_presto)), 1e-3)

})
