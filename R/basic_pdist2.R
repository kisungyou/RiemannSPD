#' Pairwise distance for two sets of data
#' 
#' 
#' @param spd1 a S3 \code{"spd"} class for \eqn{M} SPD matrices.
#' @param spd2 a S3 \code{"spd"} class for \eqn{N} SPD matrices.
#' @param geometry (case-insensitive) name of supported geometry from \code{spd.geometry("spd.pdist2")}.
#' 
#' @return an \eqn{(M\times N)} matrix of distances.
#' 
#' @concept basic
#' @export
spd.pdist <- function(spd1, spd2, geometry){
  #------------------------------------------
  # PREP
  if (!check_spdobj(spd1)){
    stop("* spd.pdist2 : 'spd1' is not a valid object of 'spd' class.")
  }
  if (!check_spdobj(spd2)){
    stop("* spd.pdist2 : 'spd2' is not a valid object of 'spd' class.")
  }
  par_geom = tolower(geometry)

  #------------------------------------------
  # COMPUTE & RETURN
  output = src_pdist2(spd1$data, spd2$data, par_geom)
}