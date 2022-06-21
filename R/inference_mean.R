#' Fréchet mean and variation
#' 
#' 
#' @param spd a S3 \code{"spd"} class for \eqn{N} of \eqn{(p\times p)} SPD matrices
#' @param geometry (case-insensitive) name of supported geometry from \code{spd.geometry("spd.mean")}.
#' @param ... extra parameters including \describe{
#' \item{weight}{a length-\eqn{N} vector of weights that sum to 1 (default: uniform).}
#' \item{maxiter}{maximum number of iterations to be run (default:100).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-8).}
#' }
#' 
#' 
#' @concept inference
#' @export
spd.mean <- function(spd, geometry, ...){
  #------------------------------------------
  # PREP
  # explicit
  if (!check_spdobj(spd)){
    stop("* spd.mean : 'spd' is not a valid object of 'spd' class.")
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
      stop("* spd.mean : 'weight' is not a valid one. Please consult the help file.")
    }
    par_weight = par_weight/base::sum(par_weight)
  } else {
    par_weight = rep(1/N, N)
  }
  
  # geometry check
  geom_avail = RiemannSPD::spd.geometry("spd.mean")
  if (!(par_geom%in%geom_avail)){
    stop(paste0("* spd.mean : the provided geometry '",par_geom,"' is not currently available."))
  }
  
  #------------------------------------------
  # COMPUTE & RETURN
  output = selector_mean(spd$data, par_weight, par_iter, par_thr, par_geom)
  return(output)
}

# # personal example
# dat3d = array(0,c(5,5,10))
# for (i in 1:10){
#   dat3d[,,i] = cov(matrix(rnorm(20*5), ncol=5))
# }
# spdobj = spd.wrap(dat3d)
# 
# spd.geometry("spd.mean")
# 
# meanAIRM = spd.mean(spdobj, "airm")
# meanLERM = spd.mean(spdobj, "lerm")
# meanEUCL = spd.mean(spdobj, "euclid")
# meanCHOL = spd.mean(spdobj, "chol")
# 
# par(mfrow=c(2,2), pty="s")
# image(meanAIRM$mean, main=paste0("AIRM:",round(meanAIRM$variation,3)))
# image(meanLERM$mean, main=paste0("LERM:",round(meanLERM$variation,3)))
# image(meanEUCL$mean, main=paste0("EUCL:",round(meanEUCL$variation,3)))
# image(meanCHOL$mean, main=paste0("CHOL:",round(meanCHOL$variation,3)))