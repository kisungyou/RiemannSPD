#' Data Harmonization by Simeon et al. (2022)
#' 
#' 
#' 
#' @references 
#' \insertRef{simeon_2022_RiemannianGeometryFunctional}{RiemannSPD}
#'  
#' @concept harmony
#' @export
spd.harmony22S <- function(spdlist, ref=c("mean","identity")){
  #------------------------------------------
  # PREP
  # explicit
  if (!check_spdlist(spdlist)){
    stop("* spd.harmony22S : 'spdlist' is not a list of valid 'spd'-class objects.")
  }
  p    = spdlist[[1]]$size[1]
  K    = length(spdlist)
  Cref = match.arg(ref)

  #------------------------------------------
  # COMPUTE
  # 1. per-class LERM as a list : each element is 3d-array
  per_class_LERM = vector(mode="list", length=K)
  for (k in 1:K){
    per_class_LERM[[k]] = aux_cube_logm(spdlist[[k]]$data)
  }
  # 2. per-class mean : (p,p,K)
  per_class_mean = array(0,c(p,p,K))
  for (k in 1:K){
    per_class_mean[,,k] = as.matrix(aux_3d_mean(per_class_LERM[[k]]))
  }
  # 3. global mean 
  if (identical(Cref,"identity")){
    global_mean = array(0,c(p,p))
  } else {
    global_mean = as.matrix(aux_3d_mean(per_class_mean))
  }
  # 4. apply
  transformed = list()
  for (k in 1:K){
    tmp_obj = spdlist[[k]]
  }
  
  # 3. apply
  transformed = list()
  for (k in 1:K){
    tmp_obj = spdlist[[k]]
    tmp_obj$data = spd.harmony22S.transform(per_class_LERM[[k]], per_class_mean[,,k], global_mean)
    transformed[[k]] = tmp_obj
  }
  
  #------------------------------------------
  # FIN
  # prepare : per-class mean
  mean_class  = vector(mode="list", length=K)
  per_class_mean_trf = aux_cube_expm(per_class_mean)
  for (k in 1:K){
    mean_class[[k]] = as.matrix(per_class_mean_trf[,,k])
  }
  # global mean computation
  mean_global = aux_expm(global_mean)
  # return
  return(list(transformed=transformed,
              mean.local=mean_class,
              mean.global=mean_global))
}

# auxiliary for harmony_22S -----------------------------------------------
#' @keywords internal
#' @noRd
spd.harmony22S.transform <- function(LE3d, classmean, globalmean){
  # prep
  p = base::nrow(globalmean)
  N = dim(LE3d)[3]
  
  # translation in T_I
  translated = array(0,c(p,p,N))
  for (n in 1:N){
    translated[,,n] = globalmean + (as.matrix(LE3d[,,n]) - classmean)
  }
  
  # return an exponentiated cube
  return(aux_cube_expm(translated))
}



# # working example
# spdlist = list()
# for (k in 1:10){
#   nsam = round(sample(20:30, 1))
#   tmper = array(0,c(5,5,nsam))
#   for (i in 1:nsam){
#     tmper[,,i] = cov(matrix(rnorm(20*5), ncol=5))
#   }
#   spdlist[[k]] = spd.wrap(tmper)
# }
# 
# #
# trflist = spd.harmony22S(spdlist)
# par(mfrow=c(3,3), pty="s")
# for (i in 1:9){
#   image(trflist$mean.class[[i]], main=paste0("class ",i))
# }
