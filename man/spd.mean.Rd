% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference_mean.R
\name{spd.mean}
\alias{spd.mean}
\title{Fréchet mean and variation}
\usage{
spd.mean(spd, geometry, ...)
}
\arguments{
\item{spd}{a S3 \code{"spd"} class for \eqn{N} of \eqn{(p\times p)} SPD matrices}

\item{geometry}{(case-insensitive) name of supported geometry from \code{spd.geometry("spd.mean")}.}

\item{...}{extra parameters including \describe{
\item{weight}{a length-\eqn{N} vector of weights that sum to 1 (default: uniform).}
\item{maxiter}{maximum number of iterations to be run (default:100).}
\item{eps}{tolerance level for stopping criterion (default: 1e-8).}
}}
}
\description{
Fréchet mean and variation
}
\concept{inference}
