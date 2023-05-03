#include <RcppArmadillo.h>
#include "casket_spd.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// AUXILIARY FUNCTIONS ---------------------------------------------------------
// [[Rcpp::export]]
arma::cube aux_cube_logm(arma::cube &data3d){
  int p = data3d.n_rows;
  int N = data3d.n_slices;
  
  arma::cube output(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    output.slice(n) = arma::real(arma::logmat_sympd(data3d.slice(n)));
  }
  return(output);
}
// [[Rcpp::export]]
arma::cube aux_cube_expm(arma::cube &data3d){
  int p = data3d.n_rows;
  int N = data3d.n_slices;
  
  arma::cube output(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    output.slice(n) = arma::real(arma::expmat_sym(data3d.slice(n)));
  }
  return(output);
}

arma::mat aux_weiszfeld(arma::cube data3d, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data3d.n_rows;
  int N = data3d.n_slices;
  
  double small_num = std::sqrt(arma::datum::eps); // check for convergence
  
  // initialize
  arma::mat y_old(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    y_old.diag() += weight(n)*data3d.slice(n).diag();
  }
  arma::mat y_new(p,p,fill::zeros);
  arma::mat y_tmp(p,p,fill::zeros);
  
  arma::vec distvec(N,fill::zeros);
  double increment = 1000.0;
  int counter = 0;
  double denominator = 0.0;
  
  // iterate
  while (increment > abstol){
    // compute the distance from a current point to data points
    for (int n=0; n<N; n++){
      distvec(n) = arma::norm(data3d.slice(n)-y_old,"fro");
    }
    if (arma::any(distvec <= small_num)){
      break;
    }
    
    // compute the denominator
    denominator = 0.0;
    for (int n=0; n<N; n++){
      denominator += weight(n)/distvec(n);
    }
    // compute the numerator
    y_tmp.fill(0.0);
    for (int n=0; n<N; n++){
      y_tmp += weight(n)*data3d.slice(n)/distvec(n);
    }
    
    // update
    y_new = y_tmp/denominator;
    increment = arma::norm(y_new-y_old,"fro");
    y_old = y_new;
    counter += 1;
    if (counter >= maxiter){
      break;
    }
  }
  return(y_old);
}
arma::mat aux_mat2mean(arma::mat x, arma::mat y){
  return(x*arma::real(arma::sqrtmat(arma::solve(x,y))));
}