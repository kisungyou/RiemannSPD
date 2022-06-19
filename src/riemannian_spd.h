#ifndef RIEMANNIAN_SPD_H_
#define RIEMANNIAN_SPD_H_

#include "RcppArmadillo.h"

// Geometry : airm (Affine-Invariant Riemannian Metric)
arma::mat  airm_exp(arma::mat x, arma::mat eta, double t);
arma::mat  airm_proj(arma::mat x, arma::mat u);
arma::mat  airm_log(arma::mat x, arma::mat y);
double     airm_dist(arma::mat x, arma::mat y);
Rcpp::List airm_mean(arma::cube data, arma::vec weight, int maxiter, double abstol);

// Geometry : lerm (Log-Euclidean Riemannian Metric)
double     lerm_dist(arma::mat x, arma::mat y);
Rcpp::List lerm_mean(arma::cube data, arma::vec weight);

// Geometry : chol (Cholesky)
double     chol_dist(arma::mat x, arma::mat y);
Rcpp::List chol_mean(arma::cube data, arma::vec weight);

// Geometry : wass (Wasserstein)

// Geometry : euclid (Euclidean)
Rcpp::List euclid_mean(arma::cube data, arma::vec weight);

// Geometry : euclid (Euclidean)
Rcpp::List euclid_mean(arma::cube data, arma::vec weight);

#endif
