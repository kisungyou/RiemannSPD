#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// ======================================================
// COLLECTION OF AUXILIARY FUNCTIONS
// aux_cube_expm    : for each slice of the cube, apply expm
// aux_cube_logm    : for each slice of the cube, apply logm
// aux_expm         : exponential map
// aux_halfinv      : S -> S^{-1/2}
// aux_3d_mean      : (p,p,N) -> (p,p) taking a mean
// aux_3d_transport : parallel transport for set of points
// ======================================================


// aux_cube_expm
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

// aux_cube_logm
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



// aux_halfinv 
// [[Rcpp::export]]
arma::mat aux_halfinv(arma::mat &X){
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, X);
  arma::mat output = eigvec*arma::diagmat(1.0/arma::sqrt(eigval))*eigvec.t();
  return(output);
}

// aux_3d_mean   : (p,p,N) -> (p,p) taking a mean
// [[Rcpp::export]]
arma::mat aux_3d_mean(arma::cube &data3d){
  arma::mat output = arma::mean(data3d,2);
  return(output);
}


// aux_3d_transport : parallel transport for set of points
// [[Rcpp::export]] 
arma::cube aux_3d_transport(arma::mat &pt_start, arma::mat &pt_end, arma::cube &data3d){
  // prep
  int p = data3d.n_rows;
  int N = data3d.n_cols;
  
  // preliminary computing
  arma::mat E = arma::sqrtmat_sympd(pt_end*arma::inv_sympd(pt_start));
  
  // iterate
  arma::cube output(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    output.slice(n) = E*data3d.slice(n)*E.t();
  }
  return(output);
}

// aux_expm         : exponential map
// [[Rcpp::export]]
arma::mat aux_expm(arma::mat &X){
  arma::mat output = arma::real(arma::expmat_sym(X));
  return(output);
}