#ifndef RCPPDIST_LST_H
#define RCPPDIST_LST_H

// NumericVector versions:

inline Rcpp::NumericVector dlst(Rcpp::NumericVector x, double df, double mu,
        double sigma, bool log_p = false){
    Rcpp::NumericVector p = Rcpp::dt((x - mu)/sigma, df, log_p);
    if ( log_p ) {
        return p - log(sigma);
    }
    else {
        return (1.0 / sigma) * p;
    }
}

inline Rcpp::NumericVector plst(Rcpp::NumericVector q, double df, double mu,
        double sigma, bool lower_tail = true, bool log_p = false){
    return Rcpp::pt((q - mu)/sigma, df, lower_tail, log_p);
}

inline Rcpp::NumericVector qlst(Rcpp::NumericVector p, double df, double mu,
        double sigma, bool lower_tail = true, bool log_p = false){
    return Rcpp::qt(p, df, lower_tail, log_p) * sigma + mu;
}

inline Rcpp::NumericVector rlst(int n, double df, double mu, double sigma){
    return Rcpp::rt(n, df) * sigma + mu;
}



// Scalar versions:

inline double d_lst(double x, double df, double mu, double sigma,
        int log_p = 0){
    double p = R::dt((x - mu) / sigma, df, log_p);
    return log_p ? (p - log(sigma)) : ((1.0 / sigma) * p);
}

inline double p_lst(double q, double df, double mu, double sigma,
        int lower_tail = 1, int log_p = 0){
    return R::pt((q - mu)/sigma, df, lower_tail, log_p);
}

inline double q_lst(double p, double df, double mu, double sigma,
        int lower_tail = 1, int log_p = 0) {
    return R::qt(p, df, lower_tail, log_p) * sigma + mu;
}

inline double r_lst(double df, double mu, double sigma){
    return R::rt(df) * sigma + mu;
}

#endif
