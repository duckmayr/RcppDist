#ifndef RCPPDIST_TRUNCLST_H
#define RCPPDIST_TRUNCLST_H

#include <lst.h>

// Vector versions:

inline Rcpp::NumericVector dtrunclst(Rcpp::NumericVector x, double df,
        double mu, double sigma, double a, double b, bool log_p = false) {
    int n = x.size();
    Rcpp::NumericVector result(n);
    double scale = p_lst(b, df, mu, sigma, 1, 0) - p_lst(a, mu, sigma, 1, 0);
    if ( log_p ) {
        scale = log(scale);
        for ( int i = 0; i < n; ++i ) {
            if ( x[i] < a || x[i] > b ) {
                result[i] = R_NegInf;
            }
            else {
                result[i] = d_lst(x[i], df, mu, sigma, 1) - scale;
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
                result[i] = d_lst(x[i], df, mu, sigma, 0) * scale;
            }
        }
    }
    return result;
}

inline Rcpp::NumericVector ptrunclst(Rcpp::NumericVector x, double df,
        double mu, double sigma, double a, double b, bool lower_tail = true,
        bool log_p = false) {
    int n = x.size();
    Rcpp::NumericVector result(n);
    double F_a = p_lst(a, df, mu, sigma, 1, 0);
    double F_b = p_lst(b, df, mu, sigma, 1, 0);
    if ( lower_tail ) {
        if ( log_p ) { // lower_tail = true AND log_p = true
            double scale = log(F_b - F_a);
            for ( int i = 0; i < n; ++i ) {
                if ( x[i] > b ) {
                    result[i] = 0.0;
                }
                else if ( x[i] < a ) {
                    result[i] = R_NegInf;
                }
                else{
                    result[i] = log(p_lst(x[i], df, mu, sigma, 1, 0) - F_a) - scale;
                }
            }
        }
        else { // lower_tail = true AND log_p = false
            double scale = 1.0 / (F_b - F_a);
            for ( int i = 0; i < n; ++i ) {
                if ( x[i] > b ) {
                    result[i] = 1.0;
                }
                else if ( x[i] < a ) {
                    result[i] = 0.0;
                }
                else {
                    result[i] = (p_lst(x[i], df, mu, sigma, 1, 0) - F_a) * scale;
                }
            }
        }
    }
    else {
        double scale = 1.0 / (F_b - F_a);
        if ( log_p ) { // lower_tail = false AND log_p = true
            for ( int i = 0; i < n; ++i ) {
                if ( x[i] > b ) {
                    result[i] = R_NegInf;
                }
                else if ( x[i] < a ) {
                    result[i] = 0.0;
                }
                else{
                    result[i] = log(1.0 - ((p_lst(x[i], df, mu, sigma, 1, 0) - F_a) * scale));
                }
            }
        }
        else { // lower_tail = false AND log_p = false
            for ( int i = 0; i < n; ++i ) {
                if ( x[i] > b ) {
                    result[i] = 0.0;
                }
                else if ( x[i] < a ) {
                    result[i] = 1.0;
                }
                else {
                    result[i] = 1.0 - ((p_lst(x[i], df, mu, sigma, 1, 0) - F_a) * scale);
                }
            }
        }
    }
    return result;
}


inline Rcpp::NumericVector qtrunclst(Rcpp::NumericVector p, double df,
        double mu, double sigma, double a, double b, bool lower_tail = true,
        bool log_p = false) {
    int n = p.size();
    Rcpp::NumericVector probs = Rcpp::clone(p);
    if ( log_p ) {
        probs = Rcpp::exp(probs);
    }
    if ( !lower_tail ) {
        probs = 1.0 - probs;
    }
    double F_a = p_lst(a, df, mu, sigma, 1, 0);
    double F_b = p_lst(b, df, mu, sigma, 1, 0);
    Rcpp::NumericVector result(n);
    for ( int i = 0; i < n; ++i ) {
        double q = q_lst(F_a + probs[i] * (F_b - F_a), df, mu, sigma, 1, 0);
        result[i] = std::min(std::max(a, q), b);
    }
    return result;
}

inline Rcpp::NumericVector rtrunclst(int n, double df, double mu, double sigma,
        double a, double b) {
    return qtrunclst(Rcpp::runif(n), df, mu, sigma, a, b);
}



// Scalar versions:

inline double d_trunclst(double x, double df, double mu, double sigma, double a,
        double b, int log_p = 0) {
    if ( x < a || x > b ) {
        return log_p ? R_NegInf : 0.0;
    }
    double scale = p_lst(b, df, mu, sigma, 1, 0) - p_lst(a, df, mu, sigma, 1, 0);
    if ( log_p ) {
        return d_lst(x, df, mu, sigma, 1) - log(scale);
    }
    return d_lst(x, df, mu, sigma, 0) / scale;
}
    
inline double p_trunclst(double x, double df, double mu, double sigma, double a,
        double b, int lower_tail = 1, int log_p = 0) {
    double F_a = p_lst(a, df, mu, sigma, 1, 0);
    double F_b = p_lst(b, df, mu, sigma, 1, 0);
    if ( lower_tail ) {
        if ( log_p ) {
            if ( x < a ) {
                return R_NegInf;
            }
            else if ( x > b ) {
                return 0.0;
            }
            else {
                return log(p_lst(x, df, mu, sigma, 1, 0) - F_a) - log(F_b - F_a);
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
                return (p_lst(x, df, mu, sigma, 1, 0) - F_a) / (F_b - F_a);
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
                return log(1.0 - ((p_lst(x, df, mu, sigma, 1, 0) - F_a) / (F_b-F_a)));
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
                return 1.0 - (p_lst(x, df, mu, sigma, 1, 0) - F_a) / (F_b - F_a);
            }
        }
    }
}

inline double q_trunclst(double p, double df, double mu, double sigma,
        double a, double b, int lower_tail = 1, int log_p = 0) {
    if ( log_p ) {
        p = exp(p);
    }
    if ( !lower_tail ) {
        p = 1.0 - p;
    }
    double F_a = p_lst(a, mu, sigma, 1, 0);
    double F_b = p_lst(b, mu, sigma, 1, 0);
    double q = q_lst(F_a + p * (F_b - F_a), df, mu, sigma, 1, 0);
    return std::min(std::max(a, q), b);
}

inline double r_trunclst(int n, double df, double mu, double sigma,
        double a, double b) {
    return q_trunclst(R::runif(0.0, 1.0), df, mu, sigma, a, b);
}

#endif
