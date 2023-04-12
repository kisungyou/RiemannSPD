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