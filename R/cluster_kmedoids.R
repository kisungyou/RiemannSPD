#' K-medoids clustering
#' 
#' Given \eqn{N} observations \eqn{\Sigma_1, \ldots, \Sigma_N} in the SPD manifold, 
#' perform \eqn{k}-medoids clustering algorithm. The computation is based on pairwise dissimilarity, which is 
#' controlled by the \code{geometry} parameter. 
#' 
#' @param spd a S3 \code{"spd"} class for \eqn{N} SPD matrices.
#' @param k the number of clusters (default: 2). If \eqn{k<2}, it returns an error.
#' @param geometry (case-insensitive) name of supported geometry from \code{spd.geometry("spd.kmedoids")}.
#' 
#' @return a named list containing\describe{
#' \item{medoids}{a length-\eqn{k} vector of medoids' indices.}
#' \item{cluster}{a length-\eqn{N} vector of class labels (from \eqn{1:k}).}
#' }
#' 
#' @concept cluster
#' @export
spd.kmedoids <- function(spd, geometry, k=2){
  #------------------------------------------
  # PREP
  par_geom = tolower(geometry)
  par_k    = max(0, round(k))
  if (par_k < 2){
    stop("* spd.kmedoids : the choice of 'k' is not valid.")
  }
  
  #------------------------------------------
  # COMPUTE
  # 1. pairwise distance
  if (inherits(spd, "dist")){
    distobj = spd
  } else {
    if (!check_spdobj(spd)){
      stop("* spd.kmedoids : 'spd' is not a valid object of 'spd' class.")
    }
    asymlist = spd.geometry.asymmetric()
    if (par_geom%in%asymlist){
      stop(paste0("* spd.kmedoids : '",geometry,"' is an asymmetric measure. Please consider using the others." ))
    }
    distobj = stats::as.dist(src_pdist(spd$data, par_geom))
  }
  
  # 3. import the function
  func.import  = utils::getFromNamespace("hidden_kmedoids", "maotai")
  obj.kmedoids = func.import(distobj, nclust=myk) 
  
  #------------------------------------------
  # FIN
  output = list()
  output$medoids = obj.kmedoids$id.med
  output$cluster = as.integer(obj.kmedoids$clustering)
  return(output)
}