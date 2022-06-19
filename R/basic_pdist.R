#' Pairwise distance
#' 
#' 
#' @param spd a S3 \code{"spd"} class for  \eqn{N} SPD matrices.
#' @param geometry (case-insensitive) name of supported geometry from \code{spd.geometry("spd.pdist")}.
#' @param as.dist logical; if \code{TRUE}, it returns \code{dist} object, else it returns a symmetric matrix. (default: \code{FALSE}).
#' 
#' @return a S3 \code{dist} object or \eqn{(N\times N)} symmetric matrix of pairwise distances according to \code{as.dist} parameter.
#' 
#' 
#' @concept basic
#' @export
spd.pdist <- function(spd, geometry, as.dist=FALSE){
  #------------------------------------------
  # PREP
  if (!check_spdobj(spd)){
    stop("* spd.pdist : 'spd' is not a valid object of 'spd' class.")
  }
  par_geom = tolower(geometry)
  par_dist = as.logical(as.dist)
  
  #------------------------------------------
  # COMPUTE
  output = src_pdist(spd$data, par_geom)
  
  #------------------------------------------
  # RETURN
  if (par_dist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}