#include <RcppArmadillo.h>

// ======================================================
// REPRESENTATION LEARNING
// (1) cpp_representation_cmds : classical MDS
// ======================================================

using namespace Rcpp;
using namespace arma;
using namespace std;



// (1) cpp_representation_cmds -------------------------------------------------
arma::mat engine_cmds(arma::mat pdist, int ndim){ // given distance matrix, return (n x ndim)
  int N = pdist.n_rows;
  arma::mat D2 = arma::pow(pdist, 2.0);
  arma::mat J  = arma::eye<arma::mat>(N,N) - (arma::ones<arma::mat>(N,N)/(static_cast<double>(N)));
  arma::mat B  = -0.5*J*D2*J;
  
  arma::vec eigval;
  arma::mat eigvec;
  
  arma::eig_sym(eigval, eigvec, B);
  arma::mat Y = eigvec.tail_cols(ndim)*arma::diagmat(arma::sqrt(eigval.tail(ndim)));
  return(Y);
}
double engine_stress(arma::mat D, arma::mat Dhat){
  int N = D.n_rows;
  
  double tobesq = 0.0;
  double term1  = 0.0; // numerator
  double term2  = 0.0; // denominator
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      tobesq = D(i,j)-Dhat(i,j);
      term1 += (tobesq*tobesq);
      term2 += D(i,j)*D(i,j);
    }
  }
  return(sqrt(term1/term2));  
}
// [[Rcpp::export]]
Rcpp::List cpp_representation_cmds(arma::mat &pdist, int ndim){
  // COMPUTE EMBEDDING & ITS DISTANCE
  int N = pdist.n_rows;
  arma::mat Y = engine_cmds(pdist, ndim);
  arma::mat DY(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      DY(i,j) = arma::norm(Y.row(i)-Y.row(j),2);
      DY(j,i) = DY(i,j);
    }
  }
  double stress = engine_stress(pdist, DY);
  
  // WRAP AND RETURN
  Rcpp::List output;
  output["embed"] = Y;
  output["stress"] = stress;
  return(output);
}