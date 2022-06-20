#' Fréchet median and variation
#' 
#' 
#' 
#' 
#' @concept inference
#' @export
spd.median <- function(spd, geometry, ...){
  #------------------------------------------
  # PREP
  # explicit
  if (!check_spdobj(spd)){
    stop("* spd.median : 'spd' is not a valid object of 'spd' class.")
  }
  par_geom = tolower(geometry)
  
  # implicit
  N      = dim(spd$data)[3] # number of observations
  params = list(...)
  pnames = names(params)
  
  if ("maxiter"%in%pnames){
    par_iter = max(5, round(params$maxiter))
  } else {
    par_iter = 50
  }
  if ("eps"%in%pnames){
    par_thr = max(.Machine$double.eps*10, params$eps)
  } else {
    par_thr = 1e-8
  }
  if ("weight"%in%pnames){
    par_weight = params$weight
    
    cond1 = (length(par_weight)==N)
    cond2 = (all(par_weight>0))
    if (cond1&&cond2){
      stop("* spd.median : 'weight' is not a valid one. Please consult the help file.")
    }
    par_weight = par_weight/base::sum(par_weight)
  } else {
    par_weight = rep(1/N, N)
  }
  
  #------------------------------------------
  # COMPUTE & RETURN
  output = selector_median(spd$data, par_weight, par_iter, par_thr, par_geom)
  return(output)
}
