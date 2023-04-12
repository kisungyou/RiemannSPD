#' Data : EEG Covariances for Event-Related Potentials
#' 
#' This dataset delivers 216 covariance matrices from EEG ERPs with 4 different 
#' known classes by types of sources. Among 60 channels, only 32 channels are 
#' taken and sample covariance matrix is computed for each participant. The 
#' data is taken from a Python library \href{https://mne.tools/stable/generated/mne.datasets.sample.data_path.html#mne.datasets.sample.data_path}{mne}'s 
#' sample data. Since the original eigenvalues are small, the data is further normalized to have unit diagonals, i.e., the correlation matrix.
#' 
#' @usage data(ERP)
#' 
#' @examples
#' \donttest{
#' ## LOAD THE DATA
#' data(ERP)
#' }
#' 
#' @format a named list containing\describe{
#' \item{spd}{an \eqn{(32\times 32\times 216)} array of SPD matrices with unit diagonals.}
#' \item{lab}{a length-\eqn{216} factor of 4 different classes.}
#' }
#' 
#' @concept data
"ERP"