// ============================================================
// OPTIMAL TRANSPORT
// (1) cpp_swdist_projection : compute the geodesic coordinates
// ============================================================

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


// (1) cpp_swdist_projection ---------------------------------------------------
arma::mat cpp_swdist_projector(int p){
  // draw : a sphere
  arma::vec theta(p,fill::randn);
  theta /= arma::norm(theta, 2);
  
  // draw : orthogonal matrix 
  arma::mat X(p,p,fill::randn);
  arma::mat Q;
  arma::mat R;
  arma::qr(Q,R,X);
  
  // return
  arma::mat output = Q*arma::diagmat(theta)*Q.t();
  return(output);
}
// [[Rcpp::export]]
arma::mat cpp_swdist_projection(arma::cube & LogX, arma::cube & LogY){
  // parameters
  int N = LogX.n_slices;
  int p = LogX.n_rows;
  arma::mat output(N,2,fill::zeros);
  
  // create a projector matrix
  arma::mat A = cpp_swdist_projector(p);
  
  // iterate
  for (int n=0; n<N; n++){
    output(n,0) = arma::trace(A*LogX.slice(n));
    output(n,1) = arma::trace(A*LogY.slice(n));
  }
  
  // return
  return(output);
}