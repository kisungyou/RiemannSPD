#' Principal Geodesic Analysis
#' 
#' @param spd a S3 \code{"spd"} class for \eqn{N} of \eqn{(p\times p)} SPD matrices
#' @param ndim an integer-valued target dimension (default: 2).
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations to be run (default:100).}
#' \item{eps}{tolerance level for stopping criterion (default: 1e-8).}
#' }
#' 
#' 
#' @concept LDR
#' @export
spd.pga <- function(spd, ndim=2, ...){
  #------------------------------------------
  # PREP
  # explicit
  par_ndim = max(1, round(ndim))
  if (!check_spdobj(spd)){
    stop("* spd.pga : 'spd' is not a valid object of 'spd' class.")
  }
  
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

  # author-defined
  par_weight = rep(1/N, N)
  par_geom = "airm"
  
  #------------------------------------------
  # COMPUTATION
  # 1. compute the Riemannian mean
  frechet_run  = selector_mean(spd$data, par_weight, par_iter, par_thr, par_geom)
  frechet_mean = frechet_run$mean
  
  # 2. tangentialize & half vectorization
  mats_tangent = aux_3d_mani2tan(frechet_mean, spd$data)
  vecs_tangent = trf_3d_vech(mats_tangent)
  
  # 3. compute mean & centering
  mean_tangent  = as.vector(base::colMeans(vecs_tangent))
  vecs_centered = array(0,c(N,ncol(vecs_tangent)))
  for (n in 1:N){
    vecs_centered[n,] = as.vector(vecs_tangent[n,])-mean_tangent
  }
  
  # 4. use the top right singular vectors V
  if (par_ndim < min(dim(vecs_centered))){
    svd_run = RSpectra::svds(vecs_centered, par_ndim)
    svd_V   = svd_run$v 
  } else if (par_ndim==min(dim(vecs_centered))){
    svd_run = base::svd(vecs_centered)
    svd_V   = svd_run$v
  } else {
    stop(paste0("* spd.pga : 'ndim' is too large. It should be <=",min(dim(vecs_centered)),"."))
  }
  
  # 5. find the low-dim coordinates
  coordinates = vecs_centered%*%svd_V
  
  # 6. reconstruction 
  recon2d   = vecs_centered%*%svd_V%*%t(svd_V)       # recon X*V*V'
  for (n in 1:N){
    recon2d[n,] = recon2d[n,] + mean_tangent         # don't forget to translate
  }
  recon3d   = trf_3d_ivech(recon2d)                  # 3d array
  recon_mat = aux_3d_tan2mani(frechet_mean, recon3d) # projection
  
  #------------------------------------------
  # RETURN
  output = list()
  output$embed = coordinates
  output$recon = recon_mat
  return(output)
}
  

