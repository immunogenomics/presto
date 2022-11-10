context('Test main presto runs on variety of input types')

library(presto)
data(y)
data(exprs)

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
    data(object_seurat)
    library(Seurat)
    object_seurat <- UpdateSeuratObject(object_seurat)
    res <- wilcoxauc(object_seurat, "cell_type")
    expect_equal(dim(res), c(40, 10))
    expect_true(all(!is.na(res)))
})


test_that('SingleCellExperiment interface works', {
    if (!requireNamespace('SingleCellExperiment', quietly = TRUE)) {
        skip('SingleCellExperiment not available')
    }

    library(SingleCellExperiment)
    data(object_sce)

    res <- wilcoxauc(object_sce, 'cell_type')
    expect_equal(dim(res), c(40, 10))
    expect_true(all(!is.na(res)))
})
