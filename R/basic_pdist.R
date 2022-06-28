#' Pairwise distance
#' 
#' Given a collection of SPD matrices \eqn{X_1, \ldots, X_N}, compute all 
#' pairwise distances under the designated geometry. 
#' 
#' @param spd a S3 \code{"spd"} class for \eqn{N} of \eqn{(p\times p)} SPD matrices.
#' @param geometry (case-insensitive) name of supported geometry from \code{spd.geometry("spd.pdist")}.
#' 
#' @return an \eqn{(N\times N)} matrix of pairwise distances.
#' 
#' @examples 
#' ## toy data
#' gen2obj = spd.gen2class(n1=5, n2=5)
#' spd_obj = gen2obj$spd     # 'spd' class object
#' 
#' ## compute distance under different geometries
#' distA = RiemannSPD::spd.pdist(spd_obj, "airm")
#' distL = spd.pdist(spd_obj, "lerm")
#' distC = spd.pdist(spd_obj, "chol")
#' distW = spd.pdist(spd_obj, "wass")
#'  
#' ## visualize different distance matrices
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(distA, main="AIRM", xaxt="n", yaxt="n")
#' image(distL, main="LERM", xaxt="n", yaxt="n")
#' image(distC, main="Cholesky", xaxt="n", yaxt="n")
#' image(distW, main="Wasserstein", xaxt="n", yaxt="n")
#' par(opar)
#' 
#' @concept basic
#' @export
spd.pdist <- function(spd, geometry){
  #------------------------------------------
  # PREP
  if (!check_spdobj(spd)){
    stop("* spd.pdist : 'spd' is not a valid object of 'spd' class.")
  }
  par_geom = tolower(geometry)

  # geometry check
  geom_avail = RiemannSPD::spd.geometry("spd.pdist")
  if (!(par_geom%in%geom_avail)){
    stop(paste0("* spd.pdist : the provided geometry '",par_geom,"' is not currently available."))
  }
  
  #------------------------------------------
  # COMPUTE
  output = src_pdist(spd$data, par_geom)
  return(output)
  
}


# # personal example
# gen2obj = spd.gen2class(n1=5, n2=5)
# spd_obj = gen2obj$spd     # 'spd' class object
# geom_list = spd.geometry("spd.pdist")
# 
# par(mfrow=c(3,4), pty="s")
# for (i in 1:10){
#   run_obj = spd.pdist(spd_obj, geom_list[i])
#   image(run_obj, main=geom_list[i], xaxt="n", yaxt="n")
# }
