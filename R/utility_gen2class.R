#' Generate 2 Types of SPD Matrices
#' 
#' This function simulates two types of covariance matrices as follows. Given a model 
#' covariance matrix \eqn{\Sigma \in \mathbf{R}^{p\times p}}, a random sample of size 
#' \eqn{3p} is generated from multivariate normal model \eqn{\mathcal{N}(0,\Sigma)} and 
#' its empirical covariance \eqn{\hat{\Sigma}} is computed as a perturbed version of the model 
#' matrix. In this function, the first model matrix is an identity matrix \eqn{\Sigma_1 = I_{p}}. 
#' The second model matrix is generated from an empirical covariance of a matrix filled with 
#' random numbers drawn from a uniform distribution \eqn{\mathcal{U}(0,1)}.
#' 
#' @param p dimensionality for SPD matrices (default: 5).
#' @param n1 the number of type 1 SPD matrices (default: 10).
#' @param n2 the number of type 2 SPD matrices (default: 10).
#' @param corr a logical; \code{TRUE} to coerce unit diagonal structure and \code{FALSE} otherwise (default: \code{FALSE}).
#' 
#' @return a named list containing\describe{
#' \item{spd}{a S3 \code{spd} class object of \eqn{(n_1+n_2)} SPD matrices.}
#' \item{lab}{a length \eqn{(n_1 + n_2)} vector of class labels of numeric values.}
#' }
#' 
#' @examples 
#' \donttest{
#' ## simple example 
#' gen2obj = spd.gen2class(n1=5, n2=5)
#' gen_spd = gen2obj$spd     # 'spd' class object
#' gen_lab = gen2obj$lab     # corresponding label
#' 
#' ## visualize all matrices
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,5), pty="s")
#' for (i in 1:10){
#'   image(gen_spd$data[,,i], xaxt="n", yaxt="n",
#'         main=paste0("class ",gen_lab[i]))
#' }
#' par(opar)
#' }
#' 
#' @concept utility
#' @export
spd.gen2class <- function(p=5, n1=10, n2=10, corr=FALSE){
  #------------------------------------------
  # PREP
  my_p  = max(2, round(p))
  my_n1 = max(1, round(n1))
  my_n2 = max(1, round(n2))
  my_nn = my_n1+my_n2
  
  par_corr = as.logical(corr)
  
  #------------------------------------------
  # COMPUTE
  # model covariances
  model_sig1 = diag(my_p)
  model_sig2 = stats::cov(matrix(stats::runif(2*my_p*my_p), ncol = my_p))
  
  # stack up
  stack_covs = array(0,c(my_p,my_p,my_nn))
  stack_labs = rep(0, my_nn)
  if (par_corr){
    for (i in 1:my_n1){
      stack_covs[,,i] = maotai::cov2corr(stats::cov(Riemann::rmvnorm(5*my_p, rep(0,my_p), model_sig1)))
      stack_labs[i] = 1
    }
    for (j in 1:my_n2){
      stack_covs[,,j+my_n1] = maotai::cov2corr(stats::cov(Riemann::rmvnorm(5*my_p, rep(0,my_p), model_sig2)))
      stack_labs[j+my_n1] = 2
    }  
  } else {
    for (i in 1:my_n1){
      stack_covs[,,i] = stats::cov(Riemann::rmvnorm(5*my_p, rep(0,my_p), model_sig1))
      stack_labs[i] = 1
    }
    for (j in 1:my_n2){
      stack_covs[,,j+my_n1] = stats::cov(Riemann::rmvnorm(5*my_p, rep(0,my_p), model_sig2))
      stack_labs[j+my_n1] = 2
    }
  }
  
  
  #------------------------------------------
  # FIN
  # wrap
  stack_spd = structure(list(data=stack_covs, size=c(my_p,my_p), corr=par_corr), class="spd")
  
  # return
  output = list()
  output$spd = stack_spd
  output$lab = stack_labs
  return(output)
}