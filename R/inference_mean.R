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
# all_geom = spd.geometry("spd.mean") # until KL, there are 9 methods.
# all_case = length(all_geom)
# 
# x11()
# par(mfrow=c(4,3), pty="s")
# for (i in 1:10){
#   image(spd.mean(spdobj, geometry=all_geom[i])$mean, xaxt='n', yaxt='n',
#         main=all_geom[i])
# }
# # compare pross with estSS from 'shapes' (works)
# dat3d = array(0,c(5,5,10))
# for (i in 1:10){
#   dat3d[,,i] = cov(matrix(rnorm(20*5), ncol=5))
# }
# spdobj = spd.wrap(dat3d)
# 
# est1 = spd.mean(spdobj, geometry="pross", eps=1e-12, maxiter=200)$mean
# est2 = shapes::estSS(dat3d)