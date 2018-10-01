// 4beta.h
//
// Copyright (C) 2018 JB Duck-Mayr
//
// This file is part of RcppDist.
//
// RcppDist is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppDist is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppDist.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RCPPDIST_4BETA_H
#define RCPPDIST_4BETA_H

// NumericVector versions:

inline Rcpp::NumericVector d4beta(const Rcpp::NumericVector& x,
        const double shape1, const double shape2, const double a,
        const double b, const bool log_p = false){
    Rcpp::NumericVector res = Rcpp::dbeta((x-a) / (b-a), shape1, shape2, log_p);
    if ( log_p ) {
        res = res - log(b - a);
        for ( int i = 0; i < x.size(); ++i ) {
            if ( x[i] < a || x[i] > b ) {
                res[i] = R_NegInf;
            }
        }
    }
    else {
        res = res / (b - a);
        for ( int i = 0; i < x.size(); ++i ) {
            if ( x[i] < a || x[i] > b ) {
                res[i] = 0.0;
            }
        }
    }
    return res;
}

inline Rcpp::NumericVector p4beta(const Rcpp::NumericVector& q,
        const double shape1, const double shape2, const double a,
        const double b, const bool lower_tail = true,
        const bool log_p = false){
    return Rcpp::pbeta((q - a) / (b - a), shape1, shape2, lower_tail, log_p);
}

inline Rcpp::NumericVector q4beta(const Rcpp::NumericVector& p,
        const double shape1, const double shape2, const double a,
        const double b, const bool lower_tail = true,
        const bool log_p = false){
    return (b - a) * Rcpp::qbeta(p, shape1, shape2, lower_tail, log_p) + a;
}

inline Rcpp::NumericVector r4beta(const int n, const double shape1,
        const double shape2, const double a, const double b){
    return (b - a) * Rcpp::rbeta(n, shape1, shape2) + a;
}



// Scalar versions:

inline double d_4beta(const double x, const double shape1, const double shape2,
        const double a, const double b, const int log_p = 0){
    if ( x < a || x > b ) {
        return log_p ? R_NegInf : 0.0;
    }
    double p = R::dbeta((x - a) / (b - a), shape1, shape2, log_p);
    return log_p ? (p - log(b - a)) : (p / (b - a));
}

inline double p_4beta(const double q, const double shape1, const double shape2,
        const double a, const double b, const int lower_tail = 1,
        const int log_p = 0){
    return R::pbeta((q - a) / (b - a), shape1, shape2, lower_tail, log_p);
}

inline double q_4beta(const double p, const double shape1, const double shape2,
        const double a, const double b, const int lower_tail = 1,
        const int log_p = 0){
    return (b - a) * R::qbeta(p, shape1, shape2, lower_tail, log_p) + a;
}

inline double r_4beta(const double shape1, const double shape2, const double a,
        const double b){
    return (b - a) * R::rbeta(shape1, shape2) + a;
}

#endif
