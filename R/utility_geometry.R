#' Show available geometries for a chosen function
#' 
#' 
#' 
#' 
#' @concept utility
#' @export
spd.geometry <- function(fname){
  if (identical(tolower(fname), "spd.pdist")){
    return(c("airm","lerm","chol","euclid"))
  } else if (identical(fname,   "spd.pdist2")){
    return(c("airm","lerm","chol","euclid"))
  } else if (identical(fname, "spd.mean")){
    return(c("airm","lerm","chol","euclid"))
  } else {
    stop("* spd.geometry : the given 'fname' does not match to one of functionalities of RiemannSPD.")
  }
}

# SUPER CLOSELY RELATED TO 'SRC_SELECTORS.CPP'.