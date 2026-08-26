context('Test main presto runs on variety of input types')

library(presto)

set.seed(42)
exprs <- matrix(rpois(25 * 150, lambda = 2), nrow = 25,
                dimnames = list(paste0("G", 1:25), NULL))
y <- rep(c("A", "B", "C"), each = 50)

test_that('presto executes on all dense and sparse 2D inputs', {
    N <- nrow(exprs) * length(unique(y))

    res_matrix <- wilcoxauc(as(exprs, 'matrix'), y)
    expect_equal(dim(res_matrix), c(N, 10))
    expect_true(all(!is.na(res_matrix)))

    res_dge <- wilcoxauc(as(exprs, 'dgeMatrix'), y)
    expect_equal(dim(res_dge), c(N, 10))
    expect_true(all(!is.na(res_dge)))

    res_dgt <- wilcoxauc(as(exprs, 'dgTMatrix'), y)
    expect_equal(dim(res_dgt), c(N, 10))
    expect_true(all(!is.na(res_dgt)))

    res_tsparse <- wilcoxauc(as(exprs, 'TsparseMatrix'), y)
    expect_equal(dim(res_tsparse), c(N, 10))
    expect_true(all(!is.na(res_tsparse)))

    res_dgc <- wilcoxauc(as(exprs, 'dgCMatrix'), y)
    expect_equal(dim(res_dgc), c(N, 10))
    expect_true(all(!is.na(res_dgc)))

    res_df <- wilcoxauc(as.data.frame(exprs), y)
    expect_equal(dim(res_df), c(N, 10))
    expect_true(all(!is.na(res_df)))
})


check_seurat <- function() {
    if (!requireNamespace('Seurat', quietly = TRUE)) {
        skip('Seurat V3 not available')
    } else {
        pkg_version <- packageVersion('Seurat')
        if (pkg_version < "3.0") {
            skip('Seurat V3 not available')
        }
    }
}

test_that("Seurat V3 interface works", {
    check_seurat()
    object_seurat <- toy_seurat()
    res <- wilcoxauc(object_seurat, "cell_type")
    expect_equal(dim(res), c(40, 10))
    expect_true(all(!is.na(res)))
})


test_that('SingleCellExperiment interface works', {
    if (!requireNamespace('SingleCellExperiment', quietly = TRUE)) {
        skip('SingleCellExperiment not available')
    }

    object_sce <- toy_sce()
    res <- wilcoxauc(object_sce, 'cell_type')
    expect_equal(dim(res), c(40, 10))
    expect_true(all(!is.na(res)))
})


test_that('toy generators are deterministic and preserve RNG state', {
    if (!requireNamespace('SingleCellExperiment', quietly = TRUE)) {
        skip('SingleCellExperiment not available')
    }

    set.seed(123)
    before <- .Random.seed
    a <- toy_sce()
    expect_identical(.Random.seed, before)

    b <- toy_sce()
    expect_identical(
        SummarizedExperiment::assay(a, 'counts'),
        SummarizedExperiment::assay(b, 'counts')
    )
})


test_that('DelayedMatrix input is processed block-wise (#26)', {
    if (!requireNamespace('DelayedArray', quietly = TRUE)) {
        skip('DelayedArray not available')
    }

    set.seed(23)
    m <- matrix(as.numeric(rpois(50 * 120, 1.5)), 50,
                dimnames = list(paste0("g", 1:50), NULL))
    m[sample(length(m), length(m) * 0.5)] <- 0
    yy <- rep(c("A", "B", "C"), 40)
    ref <- wilcoxauc(m, yy, verbose = FALSE)

    Xd <- DelayedArray::DelayedArray(m)                     # dense seed
    Xs <- DelayedArray::DelayedArray(as(m, "dgCMatrix"))    # sparse seed
    expect_equal(wilcoxauc(Xd, yy, verbose = FALSE), ref)
    expect_equal(wilcoxauc(Xs, yy, verbose = FALSE), ref)

    ## force multiple small blocks and confirm they stitch together exactly
    old <- options(presto.block.elements = 500)
    on.exit(options(old))
    expect_equal(wilcoxauc(Xd, yy, verbose = FALSE), ref)
    expect_equal(wilcoxauc(Xs, yy, verbose = FALSE), ref)

    ## NA values are caught inside a realized block
    m_na <- m
    m_na[7, 11] <- NA
    expect_error(
        wilcoxauc(DelayedArray::DelayedArray(m_na), yy, verbose = FALSE),
        "NA"
    )
})


test_that('unknown matrix classes get a helpful error', {
    expect_error(
        wilcoxauc(list(a = 1), rep("A", 3)),
        "does not know how to handle"
    )
})


test_that('SCE with a DelayedMatrix assay works end-to-end (#26)', {
    if (!requireNamespace('SingleCellExperiment', quietly = TRUE) ||
        !requireNamespace('DelayedArray', quietly = TRUE)) {
        skip('SingleCellExperiment or DelayedArray not available')
    }

    set.seed(3)
    m <- matrix(as.numeric(rpois(40 * 150, 1.2)), 40,
                dimnames = list(paste0("g", 1:40), paste0("c", 1:150)))
    m[sample(length(m), length(m) * 0.6)] <- 0
    clusters <- rep(c("A", "B", "C"), 50)

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            logcounts = DelayedArray::DelayedArray(as(m, "dgCMatrix"))
        )
    )
    sce$clusters <- clusters

    ## the exact call pattern reported in issue #26
    res <- suppressMessages(
        wilcoxauc(sce, assay = "logcounts", group_by = "clusters")
    )
    ref <- wilcoxauc(m, clusters, verbose = FALSE)
    expect_equal(res, ref)
})
