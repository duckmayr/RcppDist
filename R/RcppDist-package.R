#' RcppDist
#'
#' 'Rcpp' Integration of Additional Probability Distributions
#'
#' The Rcpp package provides a C++ library to make it easier to use C++ with R.
#' R and Rcpp provide functions for a variety of statistical distributions.
#' Several R packages make functions available to R for additional statistical distributions.
#' However, to access these functions from C++ code, a costly call to the R functions must be made.
#'
#' RcppDist provides a C++ library with functions for additional statistical distributions
#' that can be called from C++ when writing code using Rcpp or RcppArmadillo.
#' Functions are available that return NumericVectors as well as doubles,
#' and for multivariate or matrix distributions, Armadillo vectors and matrices.
#' RcppDist provides functions for the following distributions:
#' \itemize{
#'   \item The four parameter beta distribution
#'   \item The location-scale t distribution
#'   \item The truncated normal distribution
#'   \item The truncated t distribution
#'   \item A truncated location-scale t distribution
#'   \item The triangle distribution
#'   \item The multivariate normal distribution*
#'   \item The multivariate t distribution*
#'   \item The Wishart distribution*
#'   \item And the inverse Wishart distribution*.
#' }
#'
#' Distributions marked with an asterisk rely also on RcppArmadillo.
#'
#' @name RcppDist
#' @docType package
#' @author  JB Duck-Mayr
#' @useDynLib RcppDist
#' @importFrom Rcpp sourceCpp
NULL
.onUnload <- function (libpath) {
    library.dynam.unload('RcppDist', libpath)
}
