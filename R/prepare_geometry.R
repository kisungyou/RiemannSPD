#' Show available geometries for a chosen function
#' 
#' 
#' 
#' 
#' @concept prepare
#' @export
spd.geometry <- function(fname){
  if (identical(tolower(fname), "spd.pdist")){
    return(c("airm","lerm","chol","euclid","wass"))
  } else if (identical(fname,   "spd.pdist2")){
    return(c("airm","lerm","chol","euclid","wass"))
  } else if (identical(fname, "spd.mean")){
    return(c("airm","lerm","chol","euclid","wass"))
  } else if (identical(fname, "spd.median")){
    return(c("lerm","chol","euclid"))
  } else {
    stop("* spd.geometry : the given 'fname' is not provided with a variety of geometries.")
  }
}

# SUPER CLOSELY RELATED TO 'SRC_SELECTORS.CPP'.