#ifndef RIEMANNIAN_SPD_H_
#define RIEMANNIAN_SPD_H_

#include "RcppArmadillo.h"

// General Algorithm
arma::mat general_weiszfeld(arma::cube data3d, arma::vec weight, int maxiter, double abstol);

// Geometry : airm (Affine-Invariant Riemannian Metric)
arma::mat  airm_exp(arma::mat x, arma::mat eta, double t);
arma::mat  airm_proj(arma::mat x, arma::mat u);
arma::mat  airm_log(arma::mat x, arma::mat y);
double     airm_dist(arma::mat x, arma::mat y);
Rcpp::List airm_mean(arma::cube data, arma::vec weight, int maxiter, double abstol);
arma::mat  airm_transport(arma::mat x, arma::mat y, arma::mat eta);

// Geometry : lerm (Log-Euclidean Riemannian Metric)
double     lerm_dist(arma::mat x, arma::mat y);
Rcpp::List lerm_mean(arma::cube data, arma::vec weight);
Rcpp::List lerm_median(arma::cube data, arma::vec weight, int maxiter, double abstol);

// Geometry : chol (Cholesky)
double     chol_dist(arma::mat x, arma::mat y);
Rcpp::List chol_mean(arma::cube data, arma::vec weight);
Rcpp::List chol_median(arma::cube data, arma::vec weight, int maxiter, double abstol);

// Geometry : wass (2-Wasserstein)
double     wass_dist(arma::mat x, arma::mat y);
Rcpp::List wass_mean(arma::cube data, arma::vec weight, int maxiter, double abstol);

// Geometry : euclid (Euclidean)
Rcpp::List euclid_mean(arma::cube data, arma::vec weight);
Rcpp::List euclid_median(arma::cube data, arma::vec weight, int maxiter, double abstol);

// Geometry : jbld (Jensen-Bregman LogDet Divergence / S-Divergence)
double     jbld_dist(arma::mat x, arma::mat y); // square root of the divergence
Rcpp::List jbld_mean(arma::cube data, arma::vec weight, int maxiter, double abstol);
Rcpp::List jbld_median(arma::cube data, arma::vec weight, int maxiter, double abstol);

// Geometry : sqrtm (matrix square root)
double     sqrtm_dist(arma::mat x, arma::mat y);
Rcpp::List sqrtm_mean(arma::cube data, arma::vec weight);
Rcpp::List sqrtm_median(arma::cube data, arma::vec weight, int maxiter, double abstol);

// Geometry : bhat (Bhattacharyya)
double     bhat_dist(arma::mat x, arma::mat y);
Rcpp::List bhat_mean(arma::cube data, arma::vec weight, int maxiter, double abstol);
Rcpp::List bhat_median(arma::cube data, arma::vec weight, int maxiter, double abstol);

// Geometry : kl (Kullback-Leibler divergence)
double kl_dist(arma::mat x, arma::mat y);


#endif
