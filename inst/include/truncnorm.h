// truncnorm.h
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

#ifndef RCPPDIST_TRUNCNORM_H
#define RCPPDIST_TRUNCNORM_H

// Vector versions:

inline Rcpp::NumericVector dtruncnorm(const Rcpp::NumericVector& x,
        const double mu, const double sigma, const double a, const double b,
        const bool log_p = false) {
    int n = x.size();
    Rcpp::NumericVector result(n);
    double scale = R::pnorm(b, mu, sigma, 1, 0) - R::pnorm(a, mu, sigma, 1, 0);
    if ( log_p ) {
        scale = log(scale);
        for ( int i = 0; i < n; ++i ) {
            if ( x[i] < a || x[i] > b ) {
                result[i] = R_NegInf;
            }
            else {
                result[i] = R::dnorm(x[i], mu, sigma, 1) - scale;
            }
        }
    }
    else {
        scale = 1.0 / scale;
        for ( int i = 0; i < n; ++i ) {
            if ( x[i] < a || x[i] > b ) {
                result[i] = 0.0;
            }
            else {
                result[i] = R::dnorm(x[i], mu, sigma, 0) * scale;
            }
        }
    }
    return result;
}

inline Rcpp::NumericVector ptruncnorm(const Rcpp::NumericVector& x,
        const double mu, const double sigma, const double a, const double b,
        const bool lower_tail = true, const bool log_p = false) {
    int n = x.size();
    Rcpp::NumericVector result(n);
    double F_a = R::pnorm(a, mu, sigma, 1, 0);
    double F_b = R::pnorm(b, mu, sigma, 1, 0);
    if ( lower_tail ) {
        if ( log_p ) {
            double scale = log(F_b - F_a);
            for ( int i = 0; i < n; ++i ) {
                if ( x[i] > b ) {
                    result[i] = 0.0;
                }
                else if ( x[i] < a ) {
                    result[i] = R_NegInf;
                }
                else{
                result[i] = log(R::pnorm(x[i], mu, sigma, 1, 0) - F_a) - scale;
                }
            }
        }
        else {
            double scale = 1.0 / (F_b - F_a);
            for ( int i = 0; i < n; ++i ) {
                double q = std::max(std::min(x[i], b), a);
                result[i] = (R::pnorm(q, mu, sigma, 1, 0) - F_a) * scale;
            }
        }
    }
    else {
        double scale = 1.0 / (F_b - F_a);
        if ( log_p ) {
            for ( int i = 0; i < n; ++i ) {
                if ( x[i] > b ) {
                    result[i] = R_NegInf;
                }
                else if ( x[i] < a ) {
                    result[i] = 0.0;
                }
                else{
                    result[i] = log(1.0 - ((R::pnorm(x[i], mu, sigma, 1, 0) - F_a) * scale));
                }
            }
        }
        else {
            for ( int i = 0; i < n; ++i ) {
                double q = std::max(std::min(x[i], b), a);
                result[i] = 1 - ((R::pnorm(q, mu, sigma, 1, 0) - F_a) * scale);
            }
        }
    }
    return result;
}


inline Rcpp::NumericVector qtruncnorm(const Rcpp::NumericVector& p,
        const double mu, const double sigma, const double a, const double b,
        const bool lower_tail = true, const bool log_p = false) {
    int n = p.size();
    Rcpp::NumericVector probs = Rcpp::clone(p);
    if ( log_p ) {
        probs = Rcpp::exp(probs);
    }
    if ( !lower_tail ) {
        probs = 1 - probs;
    }
    double F_a = R::pnorm(a, mu, sigma, 1, 0);
    double F_b = R::pnorm(b, mu, sigma, 1, 0);
    Rcpp::NumericVector result(n);
    for ( int i = 0; i < n; ++i ) {
        double q = R::qnorm(F_a + probs[i] * (F_b - F_a), mu, sigma, 1, 0);
        result[i] = std::min(std::max(a, q), b);
    }
    return result;
}

inline Rcpp::NumericVector rtruncnorm(const int n, const double mu,
        const double sigma, const double a, const double b) {
    return qtruncnorm(Rcpp::runif(n), mu, sigma, a, b);
}



// Scalar versions:

inline double d_truncnorm(const double x, const double mu, const double sigma,
        const double a, const double b, const int log_p = 0) {
    if ( x < a || x > b ) {
        return log_p ? R_NegInf : 0.0;
    }
    double scale = R::pnorm(b, mu, sigma, 1, 0) - R::pnorm(a, mu, sigma, 1, 0);
    if ( log_p ) {
        return R::dnorm(x, mu, sigma, 1) - log(scale);
    }
    return R::dnorm(x, mu, sigma, 0) / scale;
}
    
inline double p_truncnorm(const double x, const double mu, const double sigma,
        const double a, const double b, const int lower_tail = 1,
        const int log_p = 0) {
    double F_a = R::pnorm(a, mu, sigma, 1, 0);
    double F_b = R::pnorm(b, mu, sigma, 1, 0);
    if ( lower_tail ) {
        if ( log_p ) {
            if ( x < a ) {
                return R_NegInf;
            }
            else if ( x > b ) {
                return 0.0;
            }
            else {
                return log(R::pnorm(x, mu, sigma, 1, 0) - F_a) - log(F_b - F_a);
            }
        }
        else {
            if ( x < a ) {
                return 0.0;
            }
            else if ( x > b ) {
                return 1.0;
            }
            else {
                return (R::pnorm(x, mu, sigma, 1, 0) - F_a) / (F_b - F_a);
            }
        }
    }
    else {
        if ( log_p ) {
            if ( x < a ) {
                return 0.0;
            }
            else if ( x > b ) {
                return R_NegInf;
            }
            else {
                return log(1.0 - ((R::pnorm(x, mu, sigma, 1, 0) - F_a) / (F_b-F_a)));
            }
        }
        else {
            if ( x < a ) {
                return 1.0;
            }
            else if ( x > b ) {
                return 0.0;
            }
            else {
                return 1.0 - (R::pnorm(x, mu, sigma, 1, 0) - F_a) / (F_b - F_a);
            }
        }
    }
}

inline double q_truncnorm(const double p, const double mu, const double sigma,
        const double a, const double b, const int lower_tail = 1,
        const int log_p = 0) {
    double prob = p;
    if ( log_p ) {
        prob = exp(prob);
    }
    if ( !lower_tail ) {
        prob = 1 - prob;
    }
    double F_a = R::pnorm(a, mu, sigma, 1, 0);
    double F_b = R::pnorm(b, mu, sigma, 1, 0);
    double q = R::qnorm(F_a + prob * (F_b - F_a), mu, sigma, 1, 0);
    return std::min(std::max(a, q), b);
}

inline double r_truncnorm(const double mu, const double sigma, const double a,
        const double b) {
    return q_truncnorm(R::runif(0.0, 1.0), mu, sigma, a, b);
}

#endif
