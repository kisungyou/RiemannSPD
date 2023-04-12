#' Sliced-Wasserstein Distance
#' 
#' 
#' 
#' @concept ot
#' @export
spd.swdist <- function(spd1, spd2, p=2, ...){
  # CHECK : EXPLICIT
  if (!check_spdobj(spd1)){
    stop("* spd.swdist : 'spd1' is not a valid spd-class object.")
  }
  if (!check_spdobj(spd2)){
    stop("* spd.swdist : 'spd2' is not a valid spd-class object.")
  }
  if (spd1$size[1]!=spd2$size[1]){
    stop("* spd.swdist : 'spd1' and 'spd2' are not of matching size.")
  }
  par_p   = max(1, as.double(p))
  par_dim = as.integer(spd1$size[1])
  
  # CHECK : IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("nproj"%in%pnames){
    par_nproj = max(1, round(params$nproj))
  } else {
    par_nproj = 496
  }
  
  # COMPUTE
  # log of the matrices
  LogSPD1 <- aux_cube_logm(spd1$data)
  LogSPD2 <- aux_cube_logm(spd2$data)
    
  # iterate
  rep_pdist <- rep(0, par_nproj)
  for (it in 1:par_nproj){
    # projection
    now_measures <- cpp_swdist_projection(LogSPD1, LogSPD2)
    now_measure1 <- as.vector(now_measures[,1])
    now_measure2 <- as.vector(now_measures[,2])
    
    # compute the distance
    rep_pdist[it] = R_swdist_pdist(now_measure1, now_measure2, par_p)
  }
  
  # RETURN
  return(as.double(base::mean(as.vector(rep_pdist))))
}



# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
R_swdist_pdist <- function(vec1, vec2, p){
  # convert it to ecdf
  ecdf1 = stats::ecdf(vec1)
  ecdf2 = stats::ecdf(vec2)
  
  # grid
  epsval = sqrt(.Machine$double.eps)
  seq_x  = seq(from=epsval, to=(1-epsval), length.out=1000)
  seq_y1 = as.vector(stats::quantile(ecdf1, seq_x))
  seq_y2 = as.vector(stats::quantile(ecdf2, seq_y))
  
  # compute the pdist
  ydsq   = abs(seq_y1-seq_y2)^p
  n      = length(ydsq)
  output = sum(((ydsq[1:(n-1)]+ydsq[2:n])/2)*base::diff(seq_x))
  return((output^(1/p)))
}