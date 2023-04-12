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
arma::field<arma::vec> cpp_swdist_projection(arma::cube & LogX, arma::cube & LogY){
  // parameters
  int N = LogX.n_slices;
  int M = LogY.n_slices;
  int p = LogX.n_rows;
  
  arma::vec projX(N,fill::zeros);
  arma::vec projY(M,fill::zeros);
  
  // create a projector matrix
  arma::mat A = cpp_swdist_projector(p);
  
  // iterate for X
  for (int n=0; n<N; n++){
    projX(n) = arma::trace(A*LogX.slice(n));
  }
  for (int m=0; m<M; m++){
    projY(m) = arma::trace(A*LogY.slice(m));
  }
  
  // return
  arma::field<arma::vec> output(2);
  output(0) = projX;
  output(1) = projY;
  return(output);
}