// laplace.h
//
// Copyright (C) 2019 JB Duck-Mayr
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

#ifndef RCPPDIST_LAPLACE_H
#define RCPPDIST_LAPLACE_H


// Scalar versions:

inline double d_laplace(const double x, const double mu, const double b,
                        const int log_p = 0) {
    double b_reciprocal = 1 / b;
    double scaled_abs_difference = std::abs(x - mu) * b_reciprocal;
    if ( log_p ) {
        return std::log(b_reciprocal) - M_LN2 - scaled_abs_difference;
    }
    return 0.5 * b_reciprocal * std::exp(-scaled_abs_difference);
}

inline double p_laplace(const double q, const double mu, const double b,
                        const int lower_tail = 1, const int log_p = 0) {
    double P = (q - mu) / b;
    if ( log_p ) {
        if ( lower_tail ) {
            if ( q < mu ) {
                P -= M_LN2;
            }
            else {
                P = std::log(1.0 - 0.5 * std::exp(-P));
            }
        }
        else {
            if ( q < mu ) {
                P = std::log(1.0 - 0.5 * std::exp(P));
            }
            else {
                P = -M_LN2 - P;
            }
        }
    }
    else {
        if ( lower_tail ) {
            if ( q < mu ) {
                P = 0.5 * std::exp(P);
            }
            else {
                P = 1.0 - 0.5 * std::exp(-P);
            }
        }
        else {
            if ( q < mu ) {
                P = 1.0 - 0.5 * std::exp(P);
            }
            else {
                P = 0.5 * std::exp(-P);
            }
        }
    }
    return P;
}

inline double q_laplace(const double p, const double mu, const double b,
                        const int lower_tail = 1, const int log_p = 0) {
    double P = p;
    if ( log_p ) {
        P = std::exp(P);
    }
    if ( !lower_tail ) {
        P = 1.0 - P;
    }
    P -= 0.5;
    return mu - b * std::copysign(1.0, P) * std::log(1.0 - 2.0 * std::abs(P));
}

inline double r_laplace(const double mu, const double b) {
    double U = R::runif(-0.5, 0.5);
    return mu - b * std::copysign(1.0, U) * std::log(1.0 - 2.0 * std::abs(U));
}


// NumericVector versions

inline Rcpp::NumericVector dlaplace(const Rcpp::NumericVector& x,
        const double mu, const double b, const bool log_p = false) {
    double b_reciprocal = 1 / b;
    Rcpp::NumericVector scaled_abs_difference = Rcpp::abs(x - mu) * b_reciprocal;
    if ( log_p ) {
        return std::log(b_reciprocal) - M_LN2 - scaled_abs_difference;
    }
    return 0.5 * b_reciprocal * Rcpp::exp(-scaled_abs_difference);
}

inline Rcpp::NumericVector plaplace(const Rcpp::NumericVector& q,
        const double mu, const double b, const bool lower_tail = true,
        const bool log_p = false) {
    int n = q.size();
    Rcpp::NumericVector result(n);
    for ( int i = 0; i < n; ++i ) {
        result[i] = p_laplace(q[i], mu, b, lower_tail, log_p);
    }
    return result;
}

inline Rcpp::NumericVector qlaplace(const Rcpp::NumericVector& p,
        const double mu, const double b, const bool lower_tail = true,
        const bool log_p = false) {
    int n = p.size();
    Rcpp::NumericVector result(n);
    for ( int i = 0; i < n; ++i ) {
        result[i] = q_laplace(p[i], mu, b, lower_tail, log_p);
    }
    return result;
}

inline Rcpp::NumericVector rlaplace(const int n,
        const double mu, const double b, const bool lower_tail = true,
        const bool log_p = false) {
    Rcpp::NumericVector result(n);
    for ( int i = 0; i < n; ++i ) {
        result[i] = r_laplace(mu, b);
    }
    return result;
}

#endif
