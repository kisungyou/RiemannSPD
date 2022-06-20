## AUXILIARY ROUTINES IN R
#  trf_vech     : (n,n) SPD -> n*(n+1)/2 vector : half-vectorization + diagonal
#  trf_ivech    : inverse of the operation
#  trf_3d_vech  : apply 'vech' to 3d array and stack them as rows
#  trf_3d_ivech : apply 'ivech' to row-stacked ones



# trf_vech ----------------------------------------------------------------
#' @keywords internal
#' @noRd
trf_vech <- function(symmat){
  return(c(base::diag(symmat), symmat[upper.tri(symmat)]))
}


# trf_ivech ---------------------------------------------------------------
#' @keywords internal
#' @noRd
trf_ivech <- function(symvec){
  k = length(symvec)
  p = round((-1 + sqrt(1+8*k))/2)
  
  vec1 = symvec[1:p]
  vec2 = symvec[(p+1):k]
  
  output = array(0,c(p,p))
  output[upper.tri(output)] = vec2
  output = output + t(output)
  diag(output) = vec1
  return(output)
}


# trf_3d_vech -------------------------------------------------------------
#' @keywords internal
#' @noRd
trf_3d_vech <- function(array3d){
  p = dim(array3d)[1]
  N = dim(array3d)[3]
  
  output = array(0,c(N,p*(p+1)/2))
  for (n in 1:N){
    output[n,] = trf_vech(array3d[,,n])
  }
  return(output)
}

# trf_3d_ivech ------------------------------------------------------------
#' @keywords internal
#' @noRd
trf_3d_ivech <- function(mat2d){
  N = base::nrow(mat2d)
  k = base::ncol(mat2d)
  p = round((-1 + sqrt(1+8*k))/2)
  
  output = array(0,c(p,p,N))
  for (n in 1:N){
    output[,,n] = trf_ivech(mat2d[n,])
  }
  return(output)
}