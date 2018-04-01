#ifndef RCPPDIST_RCPPDIST_H
#define RCPPDIST_RCPPDIST_H

#ifdef RCPPDIST_DONT_USE_ARMA
    #include <Rcpp.h>
#else
    #include <RcppArmadillo.h>
    #include <mvnorm.h>
    #include <mvt.h>
    #include <wishart.h>
#endif

#include <4beta.h>
#include <lst.h>
#include <truncnorm.h>
#include <trunct.h>
#include <trunclst.h>
#include <triangular.h>

#endif
