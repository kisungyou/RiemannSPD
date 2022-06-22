#include <RcppArmadillo.h>
#include "riemannian_spd.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// =============================================================================
// COLLECTION OF SELECTORS
// (1) DIST   : airm, lerm, chol, euclid, wass, jbld, sqrtm, bhat, kl
// (2) MEAN   : airm, lerm, chol, euclid, wass, jbld, sqrtm, bhat
// (3) MEDIAN :       lerm, chol, euclid,     , jbld, sqrtm, bhat
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
  } else if (geom=="jbld"){
    output = jbld_dist(x,y);
  } else if (geom=="sqrtm"){
    output = sqrtm_dist(x,y);
  } else if (geom=="bhat"){
    output = bhat_dist(x,y);
  } else if (geom=="kl"){
    output = kl_dist(x,y);
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
  } else if (geom=="jbld"){
    output = jbld_mean(data, weight, maxiter, abstol);
  } else if (geom=="sqrtm"){
    output = sqrtm_mean(data, weight);
  } else if (geom=="bhat"){
    output = bhat_mean(data, weight, maxiter, abstol);
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
  } else if (geom=="sqrtm"){
    output = sqrtm_median(data, weight, maxiter, abstol);
  } else if (geom=="jbld"){
    output = jbld_median(data, weight, maxiter, abstol);
  } else if (geom=="bhat"){
    output = bhat_median(data, weight, maxiter, abstol);
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
  arma::mat output(N,N,fill::zeros);
  if (geometry=="kl"){
    for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
        if (i!=j){
          output(i,j) = kl_dist(data.slice(i), data.slice(j));
        }
      }
    }
  } else if (geometry=="bhat"){           // special case : BHAT
    arma::vec data_det(N,fill::zeros);
    for (int n=0; n<N; n++){
      data_det(n) = arma::det(data.slice(n));
    }
    arma::mat mat_avg(p,p,fill::zeros);
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        mat_avg = (data.slice(i)+data.slice(j))/2.0;
        output(i,j) = (arma::log_det_sympd(mat_avg) - 0.5*std::log(data_det(i)*data_det(j)))*0.5;
        output(j,i) = output(i,j);
      }
    }
  } else if (geometry=="lerm"){    // special case : LERM
    arma::cube data_trf(p,p,N,fill::zeros);
    for (int n=0; n<N; n++){
      data_trf.slice(n) = arma::logmat_sympd(data.slice(n));
    }
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        output(i,j) = arma::norm(data_trf.slice(i)-data_trf.slice(j),"fro");
        output(j,i) = output(i,j);
      }
    }
  } else if (geometry=="chol"){    // special case : CHOL
    arma::cube data_trf(p,p,N,fill::zeros);
    for (int n=0; n<N; n++){
      data_trf.slice(n) = arma::chol(data.slice(n));
    }
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        output(i,j) = arma::norm(data_trf.slice(i)-data_trf.slice(j),"fro");
        output(j,i) = output(i,j);
      }
    }
  } else if (geometry=="sqrtm"){   // special case : SQRTM
    arma::cube data_trf(p,p,N,fill::zeros);
    for (int n=0; n<N; n++){
      data_trf.slice(n) = arma::sqrtmat_sympd(data.slice(n));
    }
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        output(i,j) = arma::norm(data_trf.slice(i)-data_trf.slice(j),"fro");
        output(j,i) = output(i,j);
      }
    }
  } else {                         // general cases
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        output(i,j) = selector_dist(data.slice(i), data.slice(j), geometry);
        output(j,i) = output(i,j);
      }
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
  arma::mat output(M,N,fill::zeros);
  if (geometry=="lerm"){         // special case : LERM
    arma::cube trf1(p,p,M,fill::zeros);
    arma::cube trf2(p,p,N,fill::zeros);
    
    for (int m=0; m<M; m++){
      trf1.slice(m) = arma::logmat_sympd(data1.slice(m));
    }
    for (int n=0; n<N; n++){
      trf2.slice(n) = arma::logmat_sympd(data2.slice(n));
    }
    
    for (int m=0; m<M; m++){
      for (int n=0; n<N; n++){
        output(m,n) = arma::norm(trf1.slice(m)-trf2.slice(n), "fro");
      }
    }
  } else if (geometry=="chol"){  // special case : CHOL
    arma::cube trf1(p,p,M,fill::zeros);
    arma::cube trf2(p,p,N,fill::zeros);
    
    for (int m=0; m<M; m++){
      trf1.slice(m) = arma::chol(data1.slice(m));
    }
    for (int n=0; n<N; n++){
      trf2.slice(n) = arma::chol(data2.slice(n));
    }
    
    for (int m=0; m<M; m++){
      for (int n=0; n<N; n++){
        output(m,n) = arma::norm(trf1.slice(m)-trf2.slice(n), "fro");
      }
    }
  } else if (geometry=="sqrtm"){ // special case : matrix square root
    arma::cube trf1(p,p,M,fill::zeros);
    arma::cube trf2(p,p,N,fill::zeros);
    
    for (int m=0; m<M; m++){
      trf1.slice(m) = arma::sqrtmat_sympd(data1.slice(m));
    }
    for (int n=0; n<N; n++){
      trf2.slice(n) = arma::sqrtmat_sympd(data2.slice(n));
    }
    
    for (int m=0; m<M; m++){
      for (int n=0; n<N; n++){
        output(m,n) = arma::norm(trf1.slice(m)-trf2.slice(n), "fro");
      }
    }
  } else if (geometry=="bhat"){  // special case : BHAT
    arma::vec vec_det1(M,fill::zeros);
    arma::vec vec_det2(N,fill::zeros);
    for (int m=0; m<M; m++){
      vec_det1(m) = std::real(arma::det(data1.slice(m)));
    }
    for (int n=0; n<N; n++){
      vec_det2(n) = std::real(arma::det(data2.slice(n)));
    }
    
    arma::mat mat_avg(p,p,fill::zeros);
    for (int m=0; m<M; m++){
      for (int n=0; n<N; n++){
        mat_avg     = (data1.slice(m)+data2.slice(n))/2.0;
        output(m,n) = (arma::log_det_sympd(mat_avg) - 0.5*std::log(vec_det1(m)*vec_det2(n)))*0.5;
      }
    }
  } else {                       // general cases
    for (int m=0; m<M; m++){
      for (int n=0; n<N; n++){
        output(m,n) = selector_dist(data1.slice(m), data2.slice(n), geometry);
      }
    } 
  }
  return(output);
}