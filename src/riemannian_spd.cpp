#include <RcppArmadillo.h>
#include "riemannian_spd.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


arma::mat airm_exp(arma::mat x, arma::mat eta, double t){
  arma::mat teta = t*eta;
  arma::mat tmp  = arma::expmat(arma::solve(x, teta));
  arma::mat yy   = x*tmp;
  return((yy+yy.t())/2.0);
}

arma::mat airm_proj(arma::mat x, arma::mat u){
  arma::mat outmat = (u + u.t())/2.0;
  return(outmat);
}

arma::mat airm_log(arma::mat x, arma::mat y){
  arma::mat tmp  = arma::real(arma::logmat(arma::solve(x,y)));
  arma::mat yy   = x*tmp;
  return((yy+yy.t())/2.0);
}

double airm_dist(arma::mat x, arma::mat y){
  arma::mat sol = arma::solve(x,y);
  arma::cx_mat cxXY = arma::logmat(sol);
  arma::mat XY = arma::real(cxXY);
  return(std::sqrt(arma::as_scalar(arma::trace(XY*XY))));
}

double lerm_dist(arma::mat x, arma::mat y){
  arma::mat xx = arma::real(arma::logmat_sympd(x));
  arma::mat yy = arma::real(arma::logmat_sympd(y));
  return(arma::norm(xx-yy,"fro"));
}
double chol_dist(arma::mat x, arma::mat y){
  arma::mat cx = arma::chol(x);
  arma::mat cy = arma::chol(y);
  return(arma::norm(cx-cy,"fro"));
}

Rcpp::List euclid_mean(arma::cube data, arma::vec weight){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // compute the mean
  arma::mat out_mean(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    out_mean += data.slice(n)*weight(n);
  }
  
  // compute the variation
  double dist2mean = 0.0;
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    dist2mean = arma::norm(out_mean-data.slice(n),"fro");
    out_variation += (dist2mean*dist2mean)*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["mean"] = out_mean;
  output["variation"] = out_variation;
  return(output);
}

Rcpp::List chol_mean(arma::cube data, arma::vec weight){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // compute the mean
  arma::mat tmpmat(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    tmpmat += arma::chol(data.slice(n))*weight(n);
  }
  arma::mat out_mean = arma::trans(tmpmat)*tmpmat;
  
  // compute the variation
  double dist2mean = 0.0;
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    dist2mean = chol_dist(out_mean, data.slice(n));
    out_variation += (dist2mean*dist2mean)*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["mean"] = out_mean;
  output["variation"] = out_variation;
  return(output);
}
Rcpp::List lerm_mean(arma::cube data, arma::vec weight){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // compute the mean
  arma::mat tmpmat(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    tmpmat += arma::real(arma::logmat_sympd(data.slice(n)))*weight(n);
  }
  arma::mat out_mean = arma::real(arma::expmat_sym(tmpmat));
  
  // compute the variation
  double dist2mean = 0.0;
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    dist2mean = lerm_dist(out_mean, data.slice(n));
    out_variation += (dist2mean*dist2mean)*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["mean"] = out_mean;
  output["variation"] = out_variation;
  return(output);
}
Rcpp::List airm_mean(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // stopping criterion
  double stop_inc  = 10000.0;
  int stop_counter = 0;
  
  // initialization
  arma::mat mu_old(p,p,fill::zeros);
  arma::mat mu_tmp(p,p,fill::zeros);
  arma::mat mu_new(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    mu_old.diag() += arma::diagvec(data.slice(n))*weight(n);
  }
  
  // iteration
  while (stop_inc > abstol){
    // gradient descent
    mu_tmp.fill(0.0);
    for (int n=0; n<N; n++){
      mu_tmp += weight(n)*airm_log(mu_old, data.slice(n));
    }
    mu_new = airm_exp(mu_old, mu_tmp, 1.0);
    
    // update
    stop_inc = arma::norm(mu_old-mu_new,"fro");
    mu_old   = mu_new;
    stop_counter += 1;
    if (stop_counter >= maxiter){
      break;
    }
  }
  
  // compute the variation
  double dist2mean = 0.0;
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    dist2mean = lerm_dist(mu_old, data.slice(n));
    out_variation += (dist2mean*dist2mean)*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["mean"] = mu_old;
  output["variation"] = out_variation;
  return(output);
}
