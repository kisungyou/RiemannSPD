#' RiemannSPD
#' 
#' 
#' We provide a variety of algorithms for symmetric and positive-definite (SPD) matrices, 
#' a collection of which has rich geometric characteristics.
#' 
#' @docType package
#' @name package-RiemannSPD
#' @noRd
#' @import maotai
#' @import Rdpack
#' @importFrom Riemann rmvnorm
#' @importFrom utils packageVersion getFromNamespace
#' @importFrom stats ecdf quantile rnorm cov
#' @importFrom Rcpp evalCpp
#' @useDynLib RiemannSPD, .registration=TRUE
NULL