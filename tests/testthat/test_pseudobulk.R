context('Test pseudobulk helpers')

library(presto)

test_that('collapse_counts keeps single-column meta_data a data.frame', {
    ## Regression test: meta_data[idx_keep, ] used to drop to a vector when
    ## meta_data had one column, which broke the collapse_background path of
    ## pseudobulk_deseq2().
    set.seed(1)
    m <- matrix(sample.int(5, 20 * 40, replace = TRUE), nrow = 20)
    rownames(m) <- paste0("G", 1:20)
    colnames(m) <- paste0("C", 1:40)
    meta <- data.frame(grp = sample(c("x", "y"), 40, replace = TRUE))

    cc <- collapse_counts(m, meta, "grp")
    expect_s3_class(cc$meta_data, "data.frame")
    expect_true("grp" %in% colnames(cc$meta_data))
    expect_equal(ncol(cc$counts_mat), nrow(cc$meta_data))
})

test_that('pseudobulk_deseq2 one_vs_all runs on collapsed pseudobulks', {
    if (!requireNamespace('DESeq2', quietly = TRUE)) {
        skip('DESeq2 not available')
    }
    set.seed(2)
    m <- matrix(sample.int(8, 100 * 500, replace = TRUE), nrow = 100)
    rownames(m) <- paste0("G", 1:100)
    colnames(m) <- paste0("C", 1:500)
    meta <- data.frame(
        cluster = sample(c("a", "b"), 500, replace = TRUE),
        donor = sample(paste0("d", 1:6), 500, replace = TRUE)
    )
    dc <- collapse_counts(m, meta, c("cluster", "donor"))

    res <- pseudobulk_deseq2(
        ~cluster, dc$meta_data["cluster"], dc$counts_mat,
        verbose = FALSE, present_in_min_samples = 1,
        collapse_background = FALSE, mode = "one_vs_all"
    )
    expect_s3_class(res, "data.frame")
    expect_true(all(c("group", "feature", "log2FoldChange", "padj") %in%
                    colnames(res)))
    expect_setequal(unique(res$group), c("a", "b"))
})
