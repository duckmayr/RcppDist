# RcppDist

[![Travis-CI Build Status](https://travis-ci.org/duckmayr/RcppDist.svg?branch=master)](https://travis-ci.org/duckmayr/RcppDist)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/duckmayr/RcppDist?branch=master&svg=true)](https://ci.appveyor.com/project/duckmayr/RcppDist)
[![Coverage Status](https://codecov.io/github/duckmayr/RcppDist/graph/badge.svg)](https://codecov.io/github/duckmayr/RcppDist)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

The [Rcpp package](https://github.com/RcppCore/Rcpp) provides a C++ library to make it easier to use C++ with R. R and Rcpp provide functions for a variety of statistical distributions. Several R packages make functions available to R for additional statistical distributions. However, to access these functions from C++ code, a costly call to the R functions must be made.

**RcppDist** provides a C++ library with functions for additional statistical distributions that can be called from C++ when writing code using Rcpp or [RcppArmadillo](https://github.com/RcppCore/RcppArmadillo). Functions are available that return NumericVectors as well as doubles, and for multivariate or matrix distributions, Armadillo vectors and matrices.

RcppDist will provide functions for the following distributions:
 - The four parameter beta distribution
 - The location-scale t distribution
 - The truncated normal distribution
 - The truncated t distribution
 - A truncated location-scale t distribution
 - The triangular distribution
 - The multivariate normal distribution*
 - The multivariate t distribution*
 - The Wishart distribution*
 - The inverse Wishart distribution*
 
Distributions marked with an asterisk rely on RcppArmadillo; if a user would prefer to use Rcpp but *not* RcppArmadillo (i.e. include the Rcpp headers but not the RcppArmadillo headers), include the line

    #define RCPPDIST_DONT_USE_ARMA
    
before including `RcppDist.h`, though this will make the asterisked distributions unavailable.

The distributions above were selected for inclusion because I already had occasion to generate C++ code for use with Rcpp (or RcppArmadillo) for these distributions, but I am open to requests to expand the package to include additional distributions -- just open an issue with the requested feature.

### License

GPL (>= 2)
