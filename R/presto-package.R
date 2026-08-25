#' presto: Fast differential expression
#'
#' @keywords internal
#' @useDynLib presto
#' @import Rcpp
#' @importClassesFrom Matrix dgCMatrix dgTMatrix dgeMatrix TsparseMatrix
#' @importFrom methods as is
#' @importFrom Matrix Matrix
#' @importFrom stats p.adjust pnorm rpois wilcox.test
#' @importFrom utils head
#' @importFrom Rcpp evalCpp sourceCpp loadModule
#' @importFrom rlang .data
"_PACKAGE"
