#' Data Harmonization by Yair et al. (2019)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @references 
#' \insertRef{yair_2019_harmony}{RiemannSPD}
#' 
#' @concept harmony
#' @export
spd.harmony19Y <- function(spdlist, ...){
  #------------------------------------------
  # PREP
  # explicit
  if (!check_spdlist(spdlist)){
    stop("* spd.harmony19Y : 'spdlist' is not a list of valid 'spd'-class objects.")
  }
  p = spdlist[[1]]$size[1]
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
  # 1. compute class means and stack as 3d
  class_mean_3d = array(0,c(p,p,K))
  for (k in 1:K){
    tgt_data   = spdlist[[k]]$data
    tgt_length = dim(tgt_data)[3]
    tgt_weight = rep(1.0/tgt_length, tgt_length)
    tgt_aimean = selector_mean(tgt_data, tgt_weight, par_iter, par_thr, "airm")
    
    class_mean_3d[,,k] = tgt_aimean$mean
  }
  
  # 2. compute the global mean
  global_mean = as.matrix(selector_mean(class_mean_3d, rep(1.0/K, K), par_iter, par_thr, "airm")$mean)
  
  # 3. apply parallel transport
  transformed = list()
  for (k in 1:K){
    tmp_obj = spdlist[[k]]
    tmp_obj$data = aux_3d_transport(class_mean_3d[,,k], global_mean, spdlist[[k]]$data)
    transformed[[k]] = tmp_obj
  }
  
  #------------------------------------------
  # FIN
  # prepare : per-class mean
  mean.class  = vector(mode="list", length=K)
  for (k in 1:K){
    mean.class[[k]] = as.matrix(class_mean_3d[,,k])
  }
  # return
  return(list(transformed=transformed,
              mean.local=mean.class,
              mean.global=global_mean))
}
