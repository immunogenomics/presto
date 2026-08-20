context('Test correctness of Wilcox results')
library(presto)

set.seed(42)
exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
                dimnames = list(paste0("G", 1:25), NULL))
y <- rep(c("A", "B", "C"), each = 50)

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
        .res <- res_base_r[[group]] %>% data.frame() 
        .res$feature <- row.names(.res)
        .res <- .res %>%     
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

test_that('wilcoxauc errors on NA values in X instead of silent wrong results', {
    ## Regression test for #25: NA in the data matrix used to silently
    ## corrupt the ranks and return an incorrect p-value.
    exprs_na <- exprs
    exprs_na[3, 10] <- NA

    expect_error(wilcoxauc(exprs_na, y), "NA")
    expect_error(wilcoxauc(as(exprs_na, 'dgCMatrix'), y), "NA")

    ## The clean matrix still works.
    res_ok <- wilcoxauc(exprs, y, verbose = FALSE)
    expect_equal(nrow(res_ok), nrow(exprs) * length(unique(y)))
    expect_true(all(!is.na(res_ok)))
})

test_that('sparse wilcoxauc with many zeros matches stats::wilcox.test', {
    ## Exercises the sparse rank/zero/tie handling in cpp_wilcox_stats_dgc.
    set.seed(9)
    m <- matrix(rpois(40 * 300, 1), 40, dimnames = list(paste0("s", 1:40), NULL))
    m[sample(length(m), length(m) * 0.6)] <- 0    # ~60% zeros -> many ties
    yy <- rep(c("A", "B", "C"), length.out = 300)
    N <- ncol(m)

    base <- split(seq_len(N), yy) %>% lapply(function(idx) {
        res <- Reduce(rbind, apply(m, 1, function(g) {
            wilcox.test(g[idx], g[setdiff(seq_len(N), idx)], exact = FALSE,
                        correct = TRUE) %>% broom::tidy()
        })) %>% data.frame()
        res$feature <- row.names(m)
        res
    })
    base <- Reduce(rbind, lapply(names(base), function(gr) {
        base[[gr]] %>% dplyr::mutate(group = gr)
    })) %>% dplyr::select(feature, group, statistic, p.value)

    res <- wilcoxauc(as(m, "dgCMatrix"), yy, verbose = FALSE) %>%
        dplyr::select(feature, group, statistic, p.value = pval)
    j <- dplyr::inner_join(base, res, by = c("feature", "group"),
                           suffix = c("_base", "_presto"))
    expect_lt(max(abs(j$p.value_base - j$p.value_presto)), 1e-3)
    expect_lt(max(abs(j$statistic_base - j$statistic_presto)), 1e-3)
})

test_that('dense and sparse wilcoxauc paths give identical results', {
    ## Guards the two refactored code paths (dense vs dgCMatrix) against drift.
    set.seed(11)
    m <- matrix(rpois(80 * 400, 1.2), 80, dimnames = list(paste0("g", 1:80), NULL))
    m[sample(length(m), length(m) * 0.4)] <- 0
    yy <- sample(letters[1:5], 400, replace = TRUE)

    rd <- wilcoxauc(m, yy, verbose = FALSE)
    rs <- wilcoxauc(as(m, "dgCMatrix"), yy, verbose = FALSE)
    numcols <- c("avgExpr", "logFC", "statistic", "auc",
                 "pval", "padj", "pct_in", "pct_out")
    for (cn in numcols) {
        expect_lt(max(abs(rd[[cn]] - rs[[cn]])), 1e-8)
    }
})

test_that('nthreads > 1 gives identical results to serial', {
    ## Multithreaded ranking must be a pure speedup, never a numeric change.
    set.seed(21)
    m <- matrix(rpois(60 * 500, 1), 60, dimnames = list(paste0("g", 1:60), NULL))
    m[sample(length(m), length(m) * 0.5)] <- 0
    X <- as(m, "dgCMatrix")
    yy <- sample(letters[1:4], 500, replace = TRUE)

    serial <- wilcoxauc(X, yy, verbose = FALSE, nthreads = 1)
    expect_equal(wilcoxauc(X, yy, verbose = FALSE, nthreads = 2), serial)
    expect_equal(wilcoxauc(X, yy, verbose = FALSE, nthreads = 4), serial)
})

