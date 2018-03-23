#ifndef RCPPDIST_TRUNCNORM_H
#define RCPPDIST_TRUNCNORM_H

// Vector versions:

inline Rcpp::NumericVector dtruncnorm(Rcpp::NumericVector x, double mu,
        double sigma, double a, double b, bool log_p = false) {
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

inline Rcpp::NumericVector ptruncnorm(Rcpp::NumericVector x, double mu,
        double sigma, double a, double b, bool lower_tail = true,
        bool log_p = false) {
    int n = x.size();
    Rcpp::NumericVector result(n);
    double F_a = R::pnorm(a, mu, sigma, 1, 0);
    double F_b = R::pnorm(b, mu, sigma, 1, 0);
    if ( lower_tail ) {
        if ( log_p ) {
            double scale = log(F_b - F_a);
            for ( int i = 0; i < n; ++i ) {
                double q = std::max(std::min(x[i], b), a);
                result[i] = log(R::pnorm(q, mu, sigma, 1, 0) - F_a) - scale;
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
                double q = std::max(std::min(x[i], b), a);
                result[i] = log(1-((R::pnorm(q, mu, sigma, 1, 0)-F_a) * scale));
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


inline Rcpp::NumericVector qtruncnorm(Rcpp::NumericVector p, double mu,
        double sigma, double a, double b, bool lower_tail = true,
        bool log_p = false) {
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

inline Rcpp::NumericVector rtruncnorm(int n, double mu, double sigma,
        double a, double b) {
    return qtruncnorm(Rcpp::runif(n), mu, sigma, a, b);
}



// Scalar versions:

inline double d_truncnorm(double x, double mu, double sigma, double a,
        double b, int log_p = 0) {
    if ( x < a || x > b ) {
        return log_p ? R_NegInf : 0.0;
    }
    double scale = R::pnorm(b, mu, sigma, 1, 0) - R::pnorm(a, mu, sigma, 1, 0);
    if ( log_p ) {
        return R::dnorm(x, mu, sigma, 1) - log(scale);
    }
    return R::dnorm(x, mu, sigma, 0) / scale;
}
    
inline double p_truncnorm(double x, double mu, double sigma, double a, double b,
        int lower_tail = 1, int log_p = 0) {
    double F_a = R::pnorm(a, mu, sigma, 1, 0);
    double F_b = R::pnorm(b, mu, sigma, 1, 0);
    double q = std::max(std::min(x, b), a);
    if ( lower_tail ) {
        if ( log_p ) {
            return log(R::pnorm(q, mu, sigma, 1, 0) - F_a) - log(F_b - F_a);
        }
        else {
            return (R::pnorm(q, mu, sigma, 1, 0) - F_a) / (F_b - F_a);
        }
    }
    else {
        if ( log_p ) {
            return log(1 - ((R::pnorm(q, mu, sigma, 1, 0) - F_a) / (F_b-F_a)));
        }
        else {
            return 1 - (R::pnorm(q, mu, sigma, 1, 0) - F_a) / (F_b - F_a);
        }
    }
}

inline double q_truncnorm(double p, double mu, double sigma, double a, double b,
        int lower_tail = 1, int log_p = 0) {
    if ( log_p ) {
        p = exp(p);
    }
    if ( !lower_tail ) {
        p = 1 - p;
    }
    double F_a = R::pnorm(a, mu, sigma, 1, 0);
    double F_b = R::pnorm(b, mu, sigma, 1, 0);
    double q = R::qnorm(F_a + p * (F_b - F_a), mu, sigma, 1, 0);
    return std::min(std::max(a, q), b);
}

inline double r_truncnorm(int n, double mu, double sigma, double a, double b) {
    return q_truncnorm(R::runif(0.0, 1.0), mu, sigma, a, b);
}

#endif
