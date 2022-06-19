#' Data Harmonization by Ng et al. (2014)
#' 
#' 
#' @param spdlist a length-\eqn{K} list whose elements are S3 \code{"spd"} class objects.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations to be run (default:100).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-8).}
#' }
#' 
#' @return a named list containing \describe{
#' \item{transformed}{a length-\eqn{K} list of harmonized dataset wrapped as \code{"spd"} class.}
#' \item{mean.class}{a length-\eqn{K} list of per-class/site mean.}
#' }
#' 
#' @references 
#' \insertRef{ng_2014_TransportRiemannianManifold}{RiemannSPD}
#' 
#' @concept harmony
#' @export
spd.harmony14N <- function(spdlist, ...){
  #------------------------------------------
  # PREP
  # explicit
  if (!check_spdlist(spdlist)){
    stop("* spd.harmony14N : 'spdlist' is not a list of valid 'spd'-class objects.")
  }
  K = length(spdlist)
  # implicit
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
  
  #------------------------------------------
  # COMPUTE
  # 1. compute per-class mean
  mean.class = list() 
  for (k in 1:K){
    per_nsam   = dim(spdlist[[k]]$data)[3]
    per_weight = rep(1/per_nsam, per_nsam)
    
    mean.class[[k]] = selector_mean(spdlist[[k]]$data, 
                                    per_weight, par_iter, par_thr, "airm")$mean
  }
  
  # 2. decompose
  tmp_halfinvs = list()
  for (k in 1:K){
    tmp_halfinvs[[k]] = aux_halfinv(mean.class[[k]])
  }
  
  # 3. apply
  transformed = list()
  for (k in 1:K){
    tmp_obj = spdlist[[k]]
    tmp_obj$data = whitening_array3d(spdlist[[k]]$data, tmp_halfinvs[[k]])
    transformed[[k]] = tmp_obj
  }
  
  #------------------------------------------
  # RETURN
  return(list(transformed=transformed,
              mean.class=mean.class))
}


# auxiliary ---------------------------------------------------------------
#' @keywords internal
#' @noRd
whitening_array3d <- function(data3d, multiplier){
  p = dim(data3d)[1]
  N = dim(data3d)[3]
  
  output = array(0,c(p,p,N))
  for (n in 1:N){
    output[,,n] = multiplier%*%data3d[,,n]%*%multiplier
  }
  return(output)
}