#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// ==================================================
// COLLECTION OF AUXILIARY FUNCTIONS
// aux_halfinv : S -> S^{-1/2}
// ==================================================







// aux_halfinv 
// [[Rcpp::export]]
arma::mat aux_halfinv(arma::mat &X){
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, X);
  arma::mat output = eigvec*arma::diagmat(1.0/arma::sqrt(eigval))*eigvec.t();
  return(output);
}
