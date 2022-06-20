#' Show available geometries for a chosen function
#' 
#' 
#' 
#' 
#' @concept prepare
#' @export
spd.geometry <- function(fname){
  
  # MODIFY ---------------------------------------------------------------------
  vec_dists  = c("airm","lerm","chol","euclid","wass")
  vec_mean   = c("airm","lerm","chol","euclid","wass")
  vec_median = c("lerm","chol","euclid")
  
  # ASSIGN ---------------------------------------------------------------------
  if (identical(tolower(fname), "spd.pdist")){
    return(vec_dists)
  } else if (identical(fname,   "spd.pdist2")){
    return(vec_dists)
  } else if (identical(fname, "spd.mean")){
    return(vec_mean)
  } else if (identical(fname, "spd.median")){
    return(vec_median)
  } else {
    stop("* spd.geometry : the given 'fname' is not provided with a variety of geometries.")
  }
}

# SUPER CLOSELY RELATED TO 'SRC_SELECTORS.CPP'.