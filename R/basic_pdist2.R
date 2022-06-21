#' Pairwise distance for two sets of data
#' 
#' Given two sets of SPD matrices \eqn{X_1, \ldots, X_M} and \eqn{Y_1, \ldots, Y_N}, 
#' compute all pairwise distances across two collections of observations under 
#' the designated geometry. Please see the following carefully;\describe{
#' \item{Note 1. \code{geometry="logdet"}}{
#' Since the square root of the Jensen-Bregman Log-Determinant divergence is 
#' metric, we return the square root values instead of the original divergences 
#' to be consistent with other \emph{metrics}.
#' }
#' }
#' 
#' @param spd1 a S3 \code{"spd"} class for \eqn{M} of \eqn{(p\times p)} SPD matrices.
#' @param spd2 a S3 \code{"spd"} class for \eqn{N} of \eqn{(p\times p)} SPD matrices.
#' @param geometry (case-insensitive) name of supported geometry from \code{spd.geometry("spd.pdist2")}.
#' 
#' @return an \eqn{(M\times N)} matrix of distances.
#' 
#' @examples 
#' \donttest{
#' ## toy data
#' spd_obj1 = spd.gen2class(n1=3, n2=2)$spd
#' spd_obj2 = spd.gen2class(n1=4, n2=6)$spd
#' 
#' ## compute distance under different geometries
#' distA = spd.pdist2(spd_obj1, spd_obj2, "airm")
#' distL = spd.pdist2(spd_obj1, spd_obj2, "lerm")
#' distC = spd.pdist2(spd_obj1, spd_obj2, "chol")
#' distW = spd.pdist2(spd_obj1, spd_obj2, "wass")
#'  
#' ## visualize different distance matrices
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(distA, main="AIRM", xaxt="n", yaxt="n")
#' image(distL, main="LERM", xaxt="n", yaxt="n")
#' image(distC, main="Cholesky", xaxt="n", yaxt="n")
#' image(distW, main="Wasserstein", xaxt="n", yaxt="n")
#' par(opar)
#' }
#' 
#' @concept basic
#' @export
spd.pdist2 <- function(spd1, spd2, geometry){
  #------------------------------------------
  # PREP
  if (!check_spdobj(spd1)){
    stop("* spd.pdist2 : 'spd1' is not a valid object of 'spd' class.")
  }
  if (!check_spdobj(spd2)){
    stop("* spd.pdist2 : 'spd2' is not a valid object of 'spd' class.")
  }
  par_geom = tolower(geometry)
  
  # geometry check
  geom_avail = RiemannSPD::spd.geometry("spd.pdist2")
  if (!(par_geom%in%geom_avail)){
    stop(paste0("* spd.pdist2 : the provided geometry '",par_geom,"' is not currently available."))
  }

  #------------------------------------------
  # COMPUTE & RETURN
  output = src_pdist2(spd1$data, spd2$data, par_geom)
}