#' presto
#'
#' Fast differential expression
#'
#' @name presto
#' @docType package
#' @useDynLib presto
#' @import Rcpp
#' @importClassesFrom Matrix dgCMatrix dgTMatrix dgeMatrix TsparseMatrix
#' @importFrom methods as is
#' @importFrom Matrix Matrix
#' @importFrom stats p.adjust pnorm wilcox.test
#' @importFrom utils head
#' @importFrom Rcpp evalCpp sourceCpp loadModule
#' @importFrom rlang .data
NULL
