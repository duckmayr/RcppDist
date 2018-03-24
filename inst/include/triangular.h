#ifndef RCPPDIST_TRIANGULAR_H
#define RCPPDIST_TRIANGULAR_H

// Scalar Versions

inline double d_tri(double x, double a, double b, double c, int log_p = 0) {
    if ( x < a ) {
        if ( log_p ) {
            return R_NegInf;
        }
        return 0.0;
    }
    else if ( x < c ) {
        if ( log_p ) {
            return M_LN2 + log(x - a) - log(b - a) - log(c - a);
        }
        return (2.0 * (x - a)) / ((b - a) * (c - a));
    }
    else if ( x == c ) {
        if ( log_p ) {
            return M_LN2 - log(b - a);
        }
        return 2.0 / (b - a);
    }
    else if ( x <= b ) {
        if ( log_p ) {
            return M_LN2 + log(b - x) - log(b - a) - log(b - c);
        }
        return (2.0 * (b - x)) / ((b - a) * (b - c));
    }
    else {
        if ( log_p ) {
            return R_NegInf;
        }
        return 0.0;
    }
}

inline double p_tri(double x, double a, double b, double c,
        int lower_tail = 1, int log_p = 0) {
    if ( x < a ) {
        if ( log_p ) {
            if ( lower_tail ) {
                return R_NegInf;
            }
            return 0.0;
        }
        if ( lower_tail ) {
            return 0.0;
        }
        return 1.0;
    }
    else if ( x <= c ) {
        if ( log_p ) {
            if ( lower_tail ) {
                return (2.0 * log(x - a)) - log(b - a) - log(c - a);
            }
            return log(1.0 - (pow(x-a, 2.0) / ((b-a) * (c-a))));
        }
        if ( lower_tail ) {
            return pow(x-a, 2.0) / ((b-a) * (c-a));
        }
        return 1.0 - (pow(x-a, 2.0) / ((b-a) * (c-a)));
    }
    else if ( x <= b ) {
        if ( log_p ) {
            if ( lower_tail ) {
                return log(1.0 - (pow(b-x, 2.0) / ((b-a) * (b-c))));
            }
            return (2.0 * log(b - x)) - log(b - a) - log(b - c);
        }
        if ( lower_tail ) {
            return 1.0 - (pow(b-x, 2.0) / ((b-a) * (b-c)));
        }
        return pow(b-x, 2.0) / ((b-a) * (b-c));
    }
    else {
        if ( log_p ) {
            if ( lower_tail ) {
                return 0.0;
            }
            return R_NegInf;
        }
        if ( lower_tail ) {
            return 1.0;
        }
        return 0.0;
    }
}

inline double q_tri(double p, double a, double b, double c,
        int lower_tail = 1, int log_p = 0) {
    double c_a = c - a;
    double b_a = b - a;
    if ( log_p ) {
        p = exp(p);
    }
    if ( !lower_tail ) {
        p = 1.0 - p;
    }
    if ( p < (c_a / b_a ) ) {
        return a + sqrt(b_a * c_a * p);
    }
    return b - sqrt(b_a * (b - c) * (1 - p));
}

inline double r_tri(double a, double b, double c) {
    return q_tri(R::runif(0.0, 1.0), a, b, c);
}



// NumericVector versions:

inline Rcpp::NumericVector dtri(Rcpp::NumericVector x, double a, double b,
        double c, bool log_p = false) {
    int n = x.size();
    Rcpp::NumericVector result(n);
    for ( int i = 0; i < n; ++i ) {
        result[i] = d_tri(x[i], a, b, c, log_p);
    }
    return result;
}

inline Rcpp::NumericVector ptri(Rcpp::NumericVector x, double a, double b,
        double c, bool lower_tail = true, bool log_p = false) {
    int n = x.size();
    Rcpp::NumericVector result(n);
    for ( int i = 0; i < n; ++i ) {
        result[i] = p_tri(x[i], a, b, c, lower_tail, log_p);
    }
    return result;
}

inline Rcpp::NumericVector qtri(Rcpp::NumericVector p, double a, double b,
        double c, bool lower_tail = true, bool log_p = false) {
    int n = p.size();
    Rcpp::NumericVector result(n);
    for ( int i = 0; i < n; ++i ) {
        result[i] = q_tri(p[i], a, b, c, lower_tail, log_p);
    }
    return result;
}

inline Rcpp::NumericVector rtri(int n, double a, double b, double c) {
    Rcpp::NumericVector result(n);
    for ( int i = 0; i < n; ++i ) {
        result[i] = r_tri(a, b, c);
    }
    return result;
}

#endif
