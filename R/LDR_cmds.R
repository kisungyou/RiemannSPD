#' Classical Multidimensional Scaling
#' 
#' Given \eqn{N} observations \eqn{\Sigma_1, \ldots, \Sigma_N} in the SPD manifold, 
#' apply classical multidimensional scaling to get low-dimensional representation 
#' in Euclidean space. The computation is based on pairwise dissimilarity, which is 
#' controlled by the \code{geometry} parameter. 
#' 
#' @param spd a S3 \code{"spd"} class for \eqn{N} of \eqn{(p\times p)} SPD matrices.
#' @param geometry (case-insensitive) name of supported geometry from \code{spd.geometry("spd.cmds")}.
#' @param ndim an integer-valued target dimension (default: 2).
#' 
#' @return a named list containing \describe{
#' \item{embed}{an \eqn{(N\times ndim)} matrix whose rows are embedded observations.}
#' \item{stress}{discrepancy between embedded and original distances as a measure of error.}
#' }
#' 
#' @references 
#' \insertRef{torgerson_multidimensional_1952a}{RiemannSPD}
#' 
#' @concept LDR
#' @export
spd.cmds <- function(spd, geometry, ndim=2){
  #------------------------------------------
  # PREP
  par_ndim = max(1, round(ndim))
  par_geom = tolower(geometry)
  
  #------------------------------------------
  # COMPUTE
  # 1. pairwise distance
  if (inherits(spd, "dist")){
    pdmat = as.matrix(spd)
  } else {
    if (!check_spdobj(spd)){
      stop("* spd.cmds : 'spd' is not a valid object of 'spd' class.")
    }
    pdmat = src_pdist(spd$data, par_geom)
  }
  
  # 2. classical MDS part
  output = cpp_representation_cmds(pdmat, par_ndim)
  
  #------------------------------------------
  # FIN
  return(output)
}