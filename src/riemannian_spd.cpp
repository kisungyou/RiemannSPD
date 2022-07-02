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

arma::mat general_mat2mean(arma::mat x, arma::mat y){
  return(x*arma::real(arma::sqrtmat(arma::solve(x,y))));
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
  arma::mat cx = arma::chol(x, "lower");
  arma::mat cy = arma::chol(y, "lower");
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

// square root of the divergence
double jbld_dist(arma::mat x, arma::mat y){
  std::complex<double> term2 = arma::log_det(x*y);
  double J = arma::log_det_sympd((x+y)/2.0) - 0.5*term2.real();
  return(J);
}

Rcpp::List jbld_mean(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // stopping criterion
  double stop_inc  = 10000.0;
  int stop_counter = 0;
  
  // initialization
  // Chebbi & Moaker : their arithmetic mean is enough
  arma::mat mu_old = arma::mean(data, 2);
  arma::mat mu_new(p,p,fill::zeros);
  arma::mat tmp_fixed(p,p,fill::zeros);
  
  // iterate
  while (stop_inc > abstol){
    
    tmp_fixed.fill(0.0);
    for (int n=0; n<N; n++){
      tmp_fixed += weight(n)* arma::inv_sympd((data.slice(n)+mu_old)/2.0);
    }
    mu_new = arma::inv_sympd(tmp_fixed);
    
    // update
    stop_inc = arma::norm(mu_old-mu_new,"fro");
    mu_old   = mu_new;
    
    stop_counter += 1;
    if (stop_counter >= maxiter){
      break;
    }
  }
  
  
  // compute the variation
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    out_variation += jbld_dist(mu_old, data.slice(n))*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["mean"] = mu_old;
  output["variation"] = out_variation;
  return(output);
}

double sqrtm_dist(arma::mat x, arma::mat y){
  arma::mat xx = arma::sqrtmat_sympd(x);
  arma::mat yy = arma::sqrtmat_sympd(y);
  return(arma::norm(xx-yy,"fro"));
}

Rcpp::List sqrtm_mean(arma::cube data, arma::vec weight){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // sqrtm-transform & compute
  arma::cube sqrtm_3d(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    sqrtm_3d.slice(n) = arma::sqrtmat_sympd(data.slice(n));
  }
  
  // compute mean 
  arma::mat sqrtm_mean(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    sqrtm_mean += weight(n)*sqrtm_3d.slice(n);
  }
  arma::mat output_mean = sqrtm_mean*sqrtm_mean;
  
  // compute variation
  double distval = 0.0;
  double variation = 0.0;
  for (int n=0; n<N; n++){
    distval = arma::norm(sqrtm_mean-sqrtm_3d.slice(n),"fro");
    variation += (distval*distval)*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["mean"] = output_mean;
  output["variation"] = variation;
  return(output);
}

double bhat_dist(arma::mat x, arma::mat y){
  arma::mat z = (x+y)/2.0;
  return(std::log(arma::det(z)/std::sqrt(arma::det(x)*arma::det(y)))*0.5);
}


Rcpp::List sqrtm_median(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // sqrtm transformation
  arma::cube sqrt3d(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    sqrt3d.slice(n) = arma::sqrtmat_sympd(data.slice(n));
  }
  
  // compute the median using the general algorithm
  arma::mat median_sqrtm  = general_weiszfeld(sqrt3d, weight, maxiter, abstol);
  arma::mat median_output = median_sqrtm*median_sqrtm;
  
  // compute the variation
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    out_variation += arma::norm(sqrt3d.slice(n)-median_sqrtm, "fro")*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["median"] = median_output;
  output["variation"] = out_variation;
  return(output);
}


Rcpp::List bhat_mean(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // stopping criterion
  double stop_inc  = 10000.0;
  int stop_counter = 0;
  
  // initialization
  // Chebbi & Moaker : their arithmetic mean is enough
  arma::mat mu_old = arma::mean(data, 2);
  arma::mat mu_new(p,p,fill::zeros);
  arma::mat tmp_fixed(p,p,fill::zeros);
  
  // iterate
  while (stop_inc > abstol){
    
    tmp_fixed.fill(0.0);
    for (int n=0; n<N; n++){
      tmp_fixed += weight(n)* arma::inv_sympd((data.slice(n)+mu_old)/2.0);
    }
    mu_new = arma::inv_sympd(tmp_fixed);
    
    // update
    stop_inc = arma::norm(mu_old-mu_new,"fro");
    mu_old   = mu_new;
    
    stop_counter += 1;
    if (stop_counter >= maxiter){
      break;
    }
  }
  
  // compute the variation
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    out_variation += bhat_dist(mu_old, data.slice(n))*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["mean"] = mu_old;
  output["variation"] = out_variation;
  return(output);
}

Rcpp::List jbld_median(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // stopping criterion
  double stop_inc  = 10000.0;
  int stop_counter = 0;
  double small_num = 100.0*arma::datum::eps;
  
  // initialize
  arma::mat mu_old = arma::mean(data, 2);
  arma::mat mu_new(p,p,fill::zeros);
  
  double tmp_term1 = 0.0;
  arma::mat tmp_term2(p,p,fill::zeros);
  
  arma::vec dist_to_dat(N,fill::zeros);
  
  // iterate
  while (stop_inc > abstol){
    // compute distances
    for (int n=0; n<N; n++){
      dist_to_dat(n) = std::sqrt(jbld_dist(mu_old, data.slice(n)));
    }
    if (arma::any(dist_to_dat < small_num)){
      break;
    }
    
    // intermediate quantities
    tmp_term1 = 0.0;
    for (int n=0; n<N; n++){
      tmp_term1 += weight(n)/dist_to_dat(n);
    }
    
    tmp_term2.fill(0.0);
    for (int n=0; n<N; n++){
      tmp_term2 += (weight(n)/dist_to_dat(n))*arma::inv_sympd((mu_old+data.slice(n))/2.0);
    }
    
    // update
    mu_new   = tmp_term1*arma::inv_sympd(tmp_term2);
    stop_inc = arma::norm(mu_old-mu_new,"fro");
    mu_old   = mu_new;
    
    stop_counter += 1;
    if (stop_counter >= maxiter){
      break;
    }
  }
  
  // compute the variation
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    out_variation += std::sqrt(jbld_dist(mu_old, data.slice(n)))*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["median"] = mu_old;
  output["variation"] = out_variation;
  return(output);
}

Rcpp::List bhat_median(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;
  
  // stopping criterion
  double stop_inc  = 10000.0;
  int stop_counter = 0;
  double small_num = 100.0*arma::datum::eps;
  
  // initialize
  arma::mat mu_old = arma::mean(data, 2);
  arma::mat mu_new(p,p,fill::zeros);
  
  double tmp_term1 = 0.0;
  arma::mat tmp_term2(p,p,fill::zeros);
  
  arma::vec dist_to_dat(N,fill::zeros);
  
  // iterate
  while (stop_inc > abstol){
    // compute distances
    for (int n=0; n<N; n++){
      dist_to_dat(n) = std::sqrt(jbld_dist(mu_old, data.slice(n)));
    }
    if (arma::any(dist_to_dat < small_num)){
      break;
    }
    
    // intermediate quantities
    tmp_term1 = 0.0;
    for (int n=0; n<N; n++){
      tmp_term1 += weight(n)/dist_to_dat(n);
    }
    
    tmp_term2.fill(0.0);
    for (int n=0; n<N; n++){
      tmp_term2 += (weight(n)/dist_to_dat(n))*arma::inv_sympd((mu_old+data.slice(n))/2.0);
    }
    
    // update
    mu_new   = tmp_term1*arma::inv_sympd(tmp_term2);
    stop_inc = arma::norm(mu_old-mu_new,"fro");
    mu_old   = mu_new;
    
    stop_counter += 1;
    if (stop_counter >= maxiter){
      break;
    }
  }
  
  // compute the variation
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    out_variation += std::sqrt(bhat_dist(mu_old, data.slice(n)))*weight(n);
  }
  
  // wrap & return
  Rcpp::List output;
  output["median"] = mu_old;
  output["variation"] = out_variation;
  return(output);
}

double kl_dist(arma::mat x, arma::mat y){
  return((arma::log_det_sympd(y) - arma::log_det_sympd(x) - static_cast<double>(x.n_rows) + arma::trace(arma::solve(y,x)))*0.5);
}


Rcpp::List kl_mean(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // prep
  int p = data.n_rows;
  int N = data.n_slices;

  // compute inverses & their weighted sum 
  arma::mat inv_wsum(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    inv_wsum += weight(n)*arma::inv_sympd(data.slice(n));
  }
  arma::mat kl_mean = arma::inv_sympd(inv_wsum);

  
  
  // compute the variation
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    out_variation += weight(n)*kl_dist(kl_mean, data.slice(n));
  }
  
  // wrap & return
  Rcpp::List output;
  output["mean"] = kl_mean;
  output["variation"] = out_variation;
  return(output);
}

double pross_dist(arma::mat x, arma::mat y){
  arma::mat L1  = arma::chol(x, "lower");
  arma::mat L2  = arma::chol(y, "lower");
  arma::mat L12 = L1.t()*L2;
  
  arma::mat svdW;
  arma::vec svds;
  arma::mat svdU;
  
  arma::svd(svdW, svds, svdU, L12);
  arma::mat Rhat = svdU*svdW.t();
  
  double output = arma::norm(L1 - L2*Rhat, "fro");
  return(output);
}
double pross_distchol(arma::mat L1, arma::mat L2){
  arma::mat L12 = L1.t()*L2;
  
  arma::mat svdW;
  arma::vec svds;
  arma::mat svdU;
  
  arma::svd(svdW, svds, svdU, L12);
  arma::mat Rhat = svdU*svdW.t();
  
  double output = arma::norm(L1 - L2*Rhat, "fro");
  return(output);
}


Rcpp::List pross_mean(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // params
  int p = data.n_rows;
  int N = data.n_slices;
  
  // lower-triangular cholesky
  arma::cube chol3d(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    chol3d.slice(n) = arma::chol(data.slice(n), "lower");
  }
  
  // stopping criterion
  double stop_inc  = 10000.0;
  int stop_counter = 0;
  
  // prep
  arma::mat mu_old = arma::mean(chol3d, 2);
  arma::mat mu_tmp(p,p,fill::zeros);
  arma::mat mu_new(p,p,fill::zeros);
  arma::cube stackR(p,p,N,fill::zeros);
  
  arma::mat svdW;
  arma::vec svds;
  arma::mat svdU;
  
  // iterate
  while (stop_inc > abstol){
    // update rotation matrices
    for (int n=0; n<N; n++){
      svdW.reset();
      svds.reset();
      svdU.reset();
      
      arma::svd(svdW, svds, svdU, (mu_old.t()*chol3d.slice(n)));
      stackR.slice(n) = svdU*svdW.t();
    }
    
    // update the delta
    mu_tmp.fill(0.0);
    mu_new.fill(0.0);
    for (int n=0; n<N; n++){
      mu_tmp += weight(n)*(stackR.slice(n)*chol3d.slice(n));
    }
    mu_new = arma::chol(mu_tmp*mu_tmp.t(), "lower");

    // updating information
    stop_inc = arma::norm(mu_old - mu_new, "fro");
    mu_old   = mu_new;
    
    stop_counter += 1;
    if (stop_counter >= maxiter){
      Rcpp::Rcout << "iteration " << stop_counter << " broken." << std::endl;
      break;
    }
    Rcpp::Rcout << "iteration " << stop_counter << " completed." << std::endl;
  }
  arma::mat output_mean = mu_old*mu_old.t();
  
  // variation
  // compute the variation
  double dval = 0.0;
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    dval = pross_distchol(mu_old, chol3d.slice(n));
    out_variation += weight(n)*(dval*dval);
  }
  
  // wrap & return - don't forget to take the multiplication
  Rcpp::List output;
  output["mean"] = output_mean;
  output["variation"] = out_variation;
  return(output);
}

double jeff_dist(arma::mat x, arma::mat y){
  double output = kl_dist(x,y) + kl_dist(y,x);
  return(output);
}

Rcpp::List jeff_mean(arma::cube data, arma::vec weight, int maxiter, double abstol){
  // params
  int p = data.n_rows;
  int N = data.n_slices;
  
  // compute 1 : compute an arithmetic mean
  arma::mat meanA(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    meanA += weight(n)*data.slice(n);
  }
  // compute 2 : compute a  harmonic mean
  arma::mat tmp_h(p,p,fill::zeros);
  for (int n=0; n<N; n++){
    tmp_h += weight(n)*arma::inv_sympd(data.slice(n));
  }
  arma::mat meanH = arma::inv_sympd(tmp_h);
  
  // compute : the geometric mean
  arma::mat mean_output = general_mat2mean(meanA, meanH);
  
  // compute : the variation
  double out_variation = 0.0;
  for (int n=0; n<N; n++){
    out_variation += weight(n)*jeff_dist(data.slice(n), mean_output);
  }
  
  // wrap & return - don't forget to take the multiplication
  Rcpp::List output;
  output["mean"] = mean_output;
  output["variation"] = out_variation;
  return(output);
}