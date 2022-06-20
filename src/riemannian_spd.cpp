#include <RcppArmadillo.h>
#include "riemannian_spd.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat general_weiszfeld(arma::cube data3d, arma::vec weight, int maxiter, double abstol){
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
arma::mat airm_transport(arma::mat x, arma::mat y, arma::mat eta){
  arma::mat E = arma::sqrtmat_sympd(y*arma::inv_sympd(x));
  arma::mat output = E*eta*E.t();
  return(output);
}
double wass_dist(arma::mat x, arma::mat y){
  arma::mat xhalf = arma::sqrtmat_sympd(x);
  double d2 = arma::trace(x) + arma::trace(y) - 2.0*arma::trace(arma::sqrtmat_sympd(xhalf*y*xhalf));
  return(std::sqrt(d2));
}


Rcpp::List wass_mean(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // stopping criterion
  double stop_inc  = 10000.0;
  int stop_counter = 0;
  
  // initialization
  arma::mat mu_old(p,p,fill::zeros);
  arma::mat mu_new(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    mu_old.diag() += arma::diagvec(data.slice(n))*weight(n);
  }

  // iteration
  arma::mat Skhalf(p,p,fill::zeros);
  arma::mat Skhalfinv(p,p,fill::zeros);
  arma::mat tmp_obj(p,p,fill::zeros);
  
  while (stop_inc > abstol){
    // compute the preliminary
    Skhalf.fill(0.0);
    Skhalfinv.fill(0.0);
    
    Skhalf    = arma::real(arma::sqrtmat_sympd(mu_old));
    Skhalfinv = arma::inv_sympd(Skhalf);
    
    // summation part
    tmp_obj.fill(0.0);
    for (int n=0; n<N; n++){
      tmp_obj += arma::real(arma::sqrtmat_sympd(Skhalf*data.slice(n)*Skhalf))*weight(n);
    }
    
    // compute
    mu_new = Skhalfinv*tmp_obj*Skhalfinv;
    
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
    dist2mean = wass_dist(mu_old, data.slice(n));
    out_variation += (dist2mean*dist2mean)*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["mean"] = mu_old;
  output["variation"] = out_variation;
  return(output);
}


Rcpp::List euclid_median(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int N = data.n_slices;
  
  // compute the median using the general function
  arma::mat out_median = general_weiszfeld(data, weight, maxiter, abstol);
  
  // compute the variation
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    out_variation += arma::norm(data.slice(n)-out_median,"fro")*weight(n);
  }

  // wrap & return
  Rcpp::List output;
  output["median"] = out_median;
  output["variation"] = out_variation;
  return(output);
}

Rcpp::List lerm_median(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  arma::cube datLE(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    datLE.slice(n) = arma::real(arma::logmat_sympd(data.slice(n)));
  }
  
  // compute the median using the general function
  arma::mat tmp_median = general_weiszfeld(datLE, weight, maxiter, abstol);
  arma::mat log_median = (tmp_median + tmp_median.t())/2.0;
  arma::mat out_median = arma::expmat_sym(log_median);
  
  // compute the variation
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    out_variation += arma::norm(datLE.slice(n)-log_median, "fro")*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["median"] = out_median;
  output["variation"] = out_variation;
  return(output);
}

Rcpp::List chol_median(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  arma::cube datCHOL(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    datCHOL.slice(n) = arma::chol(data.slice(n));
  }
  
  // compute the median using the general function
  arma::mat tmp_median = general_weiszfeld(datCHOL, weight, maxiter, abstol);
  arma::mat out_median = tmp_median.t()*tmp_median;
  
  // compute the variation
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    out_variation += arma::norm(datCHOL.slice(n)-tmp_median, "fro")*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["median"] = out_median;
  output["variation"] = out_variation;
  return(output);
}