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
#' @importFrom stats as.dist cov runif rnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib RiemannSPD
NULL
# pack <- "RiemannSPD"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))
