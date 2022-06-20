#include <RcppArmadillo.h>
#include "riemannian_spd.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// =============================================================================
// COLLECTION OF SELECTORS
// (1) DIST   : airm, lerm, chol, euclid, wass
// (2) MEAN   : airm, lerm, chol, euclid, wass
// (3) MEDIAN : lerm, chol, euclid
// =============================================================================

// (1) DIST --------------------------------------------------------------------
double selector_dist(arma::mat x, arma::mat y, std::string geom){
  double output = 0.0;
  if (geom=="airm"){
    output = airm_dist(x,y);
  } else if (geom=="lerm"){
    output = lerm_dist(x,y);
  } else if (geom=="chol"){
    output = chol_dist(x,y);
  } else if (geom=="euclid"){
    output = arma::norm(x-y,"fro");
  } else if (geom=="wass"){
    output = wass_dist(x,y);
  } else {
    std::string err = "* RiemannSPD : distance measure under '" + geom + "' geometry is not currently available.";
    Rcpp::stop(err);
  }
  return(output);
}

// (2) MEAN --------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List selector_mean(arma::cube &data, arma::vec &weight, int maxiter, double abstol, std::string geom){
  Rcpp::List output;
  if (geom=="airm"){
    output = airm_mean(data, weight, maxiter, abstol);
  } else if (geom=="lerm"){
    output = lerm_mean(data, weight);
  } else if (geom=="chol"){
    output = chol_mean(data, weight);
  } else if (geom=="euclid"){
    output = euclid_mean(data, weight);
  } else if (geom=="wass"){
    output = wass_mean(data, weight, maxiter, abstol);
  } else {
    std::string err = "* RiemannSPD : Frechet mean computation under '" + geom + "' geometry is not currently available.";
    Rcpp::stop(err);
  }
  return(output);
}


// (3) MEDIAN ------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List selector_median(arma::cube &data, arma::vec &weight, int maxiter, double abstol, std::string geom){
  Rcpp::List output;
  if (geom=="euclid"){
    output = euclid_median(data, weight, maxiter, abstol);
  } else if (geom=="chol"){
    output = chol_median(data, weight, maxiter, abstol);
  } else if (geom=="lerm"){
    output = lerm_median(data, weight, maxiter, abstol);
  } else {
    std::string err = "* RiemannSPD : Frechet median computation under '" + geom + "' geometry is not currently available.";
    Rcpp::stop(err);
  }
  return(output);
}










// -------- forget about the below, they are done --------------
// [[Rcpp::export]]
arma::mat src_pdist(arma::cube &data, std::string geometry){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // compute
  arma::mat mat1(p,p,fill::zeros);
  arma::mat mat2(p,p,fill::zeros);
  
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    mat1 = data.slice(i);
    for (int j=(i+1); j<N; j++){
      mat2 = data.slice(j);
      
      output(i,j) = selector_dist(mat1, mat2, geometry);
      output(j,i) = output(i,j);
    }
  }
  return(output);
}
// [[Rcpp::export]]
arma::mat src_pdist2(arma::cube &data1, arma::cube &data2, std::string geometry){
  // prep
  int p = data1.n_rows;
  int M = data1.n_slices;
  int N = data2.n_slices;
  
  // compute
  arma::mat mat1(p,p,fill::zeros);
  arma::mat mat2(p,p,fill::zeros);
  
  arma::mat output(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    mat1 = data1.slice(m);
    for (int n=0; n<N; n++){
      mat2 = data2.slice(n);
      
      output(m,n) = selector_dist(mat1, mat2, geometry);
    }
  }
  return(output);
}