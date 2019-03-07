# RcppDist

[![Travis-CI Build Status](https://travis-ci.org/duckmayr/RcppDist.svg?branch=master)](https://travis-ci.org/duckmayr/RcppDist)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/duckmayr/RcppDist?branch=master&svg=true)](https://ci.appveyor.com/project/duckmayr/RcppDist)
[![Coverage Status](https://codecov.io/github/duckmayr/RcppDist/graph/badge.svg)](https://codecov.io/github/duckmayr/RcppDist)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/RcppDist)](https://cran.r-project.org/package=RcppDist)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/last-month/RcppDist)](https://cranlogs.r-pkg.org/badges/last-month/RcppDist)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

The [Rcpp package](https://github.com/RcppCore/Rcpp) provides a C++ library to make it easier to use C++ with R. R and Rcpp provide functions for a variety of statistical distributions. Several R packages make functions available to R for additional statistical distributions. However, to access these functions from C++ code, a costly call to the R functions must be made.

**RcppDist** provides a header-only C++ library with functions for additional statistical distributions that can be called from C++ when writing code using Rcpp or [RcppArmadillo](https://github.com/RcppCore/RcppArmadillo). Functions are available that return NumericVectors as well as doubles, and for multivariate or matrix distributions, Armadillo vectors and matrices.

RcppDist provides functions for the following distributions:
 - The four parameter beta distribution
 - The Laplace distribution
 - The location-scale t distribution
 - The truncated normal distribution
 - The truncated t distribution
 - A truncated location-scale t distribution
 - The triangular distribution
 - The multivariate normal distribution
 - The multivariate t distribution
 - The Wishart distribution
 - The inverse Wishart distribution
 
## Installation

You can install RcppDist from CRAN via

```r
install.packages("RcppDist")
```

Or, you can install the development version from GitHub via

```r
remotes::install_github("duckmayr/RcppDist")
```

## Using RcppDist

### Including RcppDist headers in standalone files

You can use RcppDist in standalone C++ files using

```cpp
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
```

If you would prefer to use Rcpp but *not* RcppArmadillo (i.e. include the Rcpp headers but not the RcppArmadillo headers), instead use

```cpp
#define RCPPDIST_DONT_USE_ARMA
#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]
```

though be aware that without Armadillo, the multivariate normal, multivariate t, Wishart, and inverse Wishart distributions will be unavailable.

### Including RcppDist headers in a package

To use RcppDist in a package that links to RcppArmadillo, you must

 - Set up your package to use `RcppArmadillo`, such as via `RcppArmadillo::RcppArmadillo.package.skeleton(your_package)`
 - Add RcppDist to the LinkingTo field of your DESCRIPTION file.
 - In any C++ file that calls a `RcppDist` function, add `#include <RcppDist.h>`

To use RcppDist in a package that does not link to RcppArmadillo, you must

 - Set up your package to use `Rcpp`, such as via `Rcpp::Rcpp.package.skeleton(your_package)` or `devtools::use_rcpp(your_package)`.
 - Add RcppDist to the LinkingTo field of your DESCRIPTION file.
 - In any C++ file that calls a `RcppDist` function, add `#include <RcppDist.h>`.
 - Use `#define RCPPDIST_DONT_USE_ARMA` before any include of `RcppDist.h`.

### RcppDist functions

Much like distributions in R, functions are prefixed by d, p, q, and r to mean density, distribution, quantile, and random number generating functions respectively. Functions that return a double rather than, say, a NumericVector are instead prefixed by d_, p_, q_, and r_. For example,

```cpp
d4beta(x, 2.0, 2.0, -5.0, 5.0)
```

gives the density of the four-parameter beta distribution with shape values 2 and 2 defined on the interval [-5, 5] at the values in the `Rcpp::NumericVector` x (the function's return value is also a `Rcpp::NumericVector`), while

```
d_4beta(x, 2.0, 2.0, -5.0, 5.0)
```

gives the density of the four-parameter beta distribution with shape values 2 and 2 defined on the interval [-5, 5] at the value in the `double` x (the function's return value is also a `double`).

Definitions and descriptions of the C++ functions provided by RcppDist, as well as more information on including RcppDist headers (such as using headers for only one or more specific distributions), are given in the vignette, which you can access from R using

```r
vignette("RcppDist")
```

An example of using RcppDist's multivariate normal generator can be seen in the `bayeslm()` R-facing function; its code is displayed in the help file for it, accessible via

```r
help("bayeslm")
```

## Contributing and requests

The distributions above were selected for inclusion because I already had occasion to generate C++ code for use with Rcpp (or RcppArmadillo) for these distributions, but I am open to requests to expand the package to include additional distributions -- just open an issue with the requested feature, or feel free to contribute the code yourself and open a pull request.

## License

GPL (>= 2)

