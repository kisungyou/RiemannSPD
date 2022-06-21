#' Show available geometries
#' 
#' @description this is description.
#' 
#' @param fname a function name for which the list of available geometries should be presented.
#' @return a vector of geometry names that are available for the given function.
#' 
#' @section Available Geometries:
#' \describe{
#' \item{\code{"airm"}}{(\emph{Affine-Invariant Riemannian Metric}) - .}
#' \item{\code{"bhat"}}{Bhattacharyya distance. not a metric.}
#' \item{\code{"chol"}}{Cholesky. See \insertCite{wang_2004_ConstrainedVariationalPrinciple}{RiemannSPD}.}
#' \item{\code{"euclid"}}{Euclidean}
#' \item{\code{"jbld"}}{(\emph{Jensen-Bregman Log-Determinant divergence}) -  
#' jensen-bregman log-determinant devivergence. also called as S-divergence from 
#' Sra. when used for distance computation, its square root is returned. Mean 
#' computation by Chebbi & Moakher. Originally by Cherian 2011. See burgeoning-033.
#' }
#' \item{\code{"lerm"}}{(\emph{Log-Euclidean Riemannian Metric}) - }
#' \item{\code{"sqrtm"}}{(\emph{Matrix Square Root}) - It was proposed in 
#' \insertCite{dryden_2009_NonEuclideanStatisticsCovariance;textual}{RiemannSPD} without 
#' much discussion at rigor. This geometry considers square root of matrices 
#' and measures distances of those under the Frobenius norm.
#' }
#' \item{\code{"wass"}}{2-wasserstein.}
#' }
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept prepare
#' @export
spd.geometry <- function(fname){
  # MODIFY ---------------------------------------------------------------------
  vec_dists  = c("airm","lerm","chol","euclid","wass","jbld","sqrtm","bhat")
  vec_mean   = c("airm","lerm","chol","euclid","wass","jbld","sqrtm")
  vec_median = c("lerm","chol","euclid","sqrtm")
  
  # vec_all = union(union(vec_dists, vec_mean), vec_median)
  # print(sort(vec_all))
  
  # ASSIGN ---------------------------------------------------------------------
  if (identical(tolower(fname), "spd.pdist")){
    return(vec_dists)
  } else if (identical(fname, "spd.pdist2")){
    return(vec_dists)
  } else if (identical(fname, "spd.mean")){
    return(vec_mean)
  } else if (identical(fname, "spd.median")){
    return(vec_median)
  } else if (identical(fname, "spd.cmds")){
    return(vec_dists)
  } else if (identical(fname, "spd.kmedoids")){
    return(vec_dists)
  } else {
    stop("* spd.geometry : the given 'fname' is not provided with a variety of geometries.")
  }
}

# SUPER CLOSELY RELATED TO 'SRC_SELECTORS.CPP'.