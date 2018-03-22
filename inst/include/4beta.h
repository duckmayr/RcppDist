#ifndef RCPPDIST_4BETA_H
#define RCPPDIST_4BETA_H

// NumericVector versions:

inline Rcpp::NumericVector d4beta(Rcpp::NumericVector x, double shape1,
        double shape2, double a, double b, bool log_p = false){
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

inline Rcpp::NumericVector p4beta(Rcpp::NumericVector q, double shape1,
        double shape2, double a, double b, bool lower_tail = true,
        bool log_p = false){
    return Rcpp::pbeta((q - a) / (b - a), shape1, shape2, lower_tail, log_p);
}

inline Rcpp::NumericVector q4beta(Rcpp::NumericVector p, double shape1,
        double shape2, double a, double b, bool lower_tail = true,
        bool log_p = false){
    return (b - a) * Rcpp::qbeta(p, shape1, shape2, lower_tail, log_p) + a;
}

inline Rcpp::NumericVector r4beta(int n, double shape1, double shape2,
        double a, double b){
    return (b - a) * Rcpp::rbeta(n, shape1, shape2) + a;
}



// Scalar versions:

inline double d_4beta(double x, double shape1, double shape2, double a,
        double b, int log_p = 0){
    if ( x < a || x > b ) {
        return log_p ? R_NegInf : 0.0;
    }
    double p = R::dbeta((x - a) / (b - a), shape1, shape2, log_p);
    return log_p ? (p - log(b - a)) : (p / (b - a));
}

inline double p_4beta(double q, double shape1, double shape2,
        double a, double b, int lower_tail = 1, int log_p = 0){
    return R::pbeta((q - a) / (b - a), shape1, shape2, lower_tail, log_p);
}

inline double q_4beta(double p, double shape1, double shape2,
        double a, double b, int lower_tail = 1, int log_p = 0){
    return (b - a) * R::qbeta(p, shape1, shape2, lower_tail, log_p) + a;
}

inline double r_4beta(double shape1, double shape2, double a, double b){
    return (b - a) * R::rbeta(shape1, shape2) + a;
}

#endif
