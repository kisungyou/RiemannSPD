#' Prepare Data on the SPD Manifold
#' 
#' The collection of symmetric positive-definite matrices is a well-known example 
#' of matrix manifold. It is defined as
#' \deqn{\mathcal{S}_{++}^p = \lbrace X \in \mathbf{R}^{p\times p} ~\vert~ X^\top = X,~ \textrm{rank}(X)=p \rbrace}
#' where the rank condition means it is strictly positive definite. In the \pkg{RiemannSPD}, 
#' we have a very strict policy to work with PD matrices only.
#' 
#' @param input SPD data matrices to be wrapped as \code{spd} class. Following inputs are considered,
#' \describe{
#' \item{array}{a \eqn{(p\times p\times n)} array where each slice along 3rd dimension is a SPD matrix.}
#' \item{list}{a length-\eqn{n} list whose elements are \eqn{(p\times p)} SPD matrices.}
#' }
#' @param corr a logical; \code{TRUE} to coerce unit diagonal structure and \code{FALSE} otherwise (default: \code{FALSE}).
#' 
#' @return a named \code{spd} S3 object containing
#' \describe{
#'   \item{data}{a \eqn{(p\times p\times n)} 3d array of SPD matrices.}
#'   \item{size}{size of each SPD matrix.}
#'   \item{corr}{a logical for correlation-valued objects.}
#' }
#' 
#' @examples
#' \donttest{
#' # LOAD THE DATA
#' data(ERP)
#' 
#' # WRAP THE DATA
#' spdobj <- spd.wrap(ERP$spd)
#' }
#' 
#' @concept utility
#' @export
spd.wrap <- function(input, corr=FALSE){
  # take either 3d array of a list
  if (is.array(input)){
    if (!check_3darray(input, symmcheck=TRUE)){
      stop("* spd.wrap : 'input' does not follow the size requirement as described.")
    }
    N = dim(input)[3]
    tmpdata = list()
    for (n in 1:N){
      tmpdata[[n]] = input[,,n]
    }
  } else if (is.list(input)){
    tmpdata = input
  } else {
    stop("* spd.wrap : 'input' should be either a 3d array or a list.")
  }
  
  # check all same size
  if (!check_list_eqsize(tmpdata, check.square=TRUE)){
    stop("* spd.wrap : elements are not of same size.")
  }
  # check a valid spd object
  N = length(tmpdata)
  for (n in 1:N){
    tmpdata[[n]] = check_spd(tmpdata[[n]], n)
  } 
  # coerce to the correlation
  par_corr = as.logical(corr)
  
  # return
  output = list()
  output$data = spdlist2array(tmpdata, par_corr)
  output$size = dim(tmpdata[[1]])
  output$corr = par_corr
  return(structure(output, class="spd"))  
}


# checkers ----------------------------------------------------------------
#' @keywords internal
#' @noRd
check_3darray <- function(x, symmcheck=TRUE){
  cond1 = is.array(x)
  cond2 = (length(dim(x))==3)
  if (symmcheck){
    cond3 = (dim(x)[1] == dim(x)[2])  
  } else {
    cond3 = TRUE
  }
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
#' @keywords internal
#' @noRd
check_list_eqsize <- function(dlist, check.square=FALSE){
  if (is.vector(dlist[[1]])){
    cond0 = all(unlist(lapply(dlist, is.vector))==TRUE)        # all vectors
    cond1 = (length(unique(unlist(lapply(dlist, length))))==1) # same length
    if (cond0&&cond1){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    cond0 = all(unlist(lapply(dlist, is.matrix))==TRUE)      # all matrices
    cond1 = (length(unique(unlist(lapply(dlist, nrow))))==1) # same row size
    cond2 = (length(unique(unlist(lapply(dlist, ncol))))==1) # same col size
    if (check.square){
      cond3 = (nrow(dlist[[1]])==ncol(dlist[[1]]))
    } else {
      cond3 = TRUE
    }
    if (cond0&&cond1&&cond2&&cond3){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}
#' @keywords internal
#' @noRd
check_spd <- function(x, id){
  p = nrow(x)
  cond1 = (nrow(x)==ncol(x))
  cond2 = (min(base::eigen(x, only.values=TRUE)$values)>.Machine$double.eps)
  cond3 = isSymmetric(x)
  if (cond1&&cond2&&cond3){
    return(x)
  } else {
    stop(paste0(" spd.wrap : the ",id,"st/nd/th object is not a valid SPD object."))
  }
}
#' @keywords internal
#' @noRd
spdlist2array <- function(datlist, coerce_corr){
  p = base::nrow(datlist[[1]])
  N = length(datlist)
  
  output = array(0,c(p,p,N))
  if (coerce_corr){
    for (n in 1:N){
      output[,,n] = maotai::cov2corr(datlist[[n]])
    }
  } else {
    for (n in 1:N){
      output[,,n] = datlist[[n]]
    } 
  }
  return(output)
}
#' @keywords internal
#' @noRd
check_spdobj <- function(spdobj){
  if (!inherits(spdobj, "spd")){
    return(FALSE)
  }
  if (!is.list(spdobj)){
    return(FALSE)
  }
  fields = names(spdobj)
  if (!("data"%in%fields)){
    return(FALSE)
  }
  if (!("size"%in%fields)){
    return(FALSE)
  }
  if (!("corr"%in%fields)){
    return(FALSE)
  }
  return(TRUE)
}
#' @keywords internal
#' @noRd
check_spdlist <- function(spdlist){
  if (!is.list(spdlist)){
    return(FALSE)
  }
  K = length(spdlist)
  for (k in 1:K){
    if (!check_spdobj(spdlist[[k]])){
      return(FALSE)
    }
  }
  vec_dim = rep(0,K)
  for (k in 1:K){
    vec_dim[k] = dim(spdlist[[k]]$data)[1]
  }
  if (!(length(unique(vec_dim))==1)){
    return(FALSE)
  }
  return(TRUE)
}
