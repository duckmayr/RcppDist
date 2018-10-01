#ifndef RCPPDIST_LST_H
#define RCPPDIST_LST_H

// NumericVector versions:

inline Rcpp::NumericVector dlst(const Rcpp::NumericVector& x, const double df,
        const double mu, const double sigma, const bool log_p = false){
    Rcpp::NumericVector p = Rcpp::dt((x - mu)/sigma, df, log_p);
    if ( log_p ) {
        return p - log(sigma);
    }
    else {
        return (1.0 / sigma) * p;
    }
}

inline Rcpp::NumericVector plst(const Rcpp::NumericVector& q, const double df,
        const double mu, const double sigma, const bool lower_tail = true,
        const bool log_p = false){
    return Rcpp::pt((q - mu)/sigma, df, lower_tail, log_p);
}

inline Rcpp::NumericVector qlst(const Rcpp::NumericVector& p, const double df,
        const double mu, const double sigma, const bool lower_tail = true,
        const bool log_p = false){
    return Rcpp::qt(p, df, lower_tail, log_p) * sigma + mu;
}

inline Rcpp::NumericVector rlst(const int n, const double df, const double mu,
        const double sigma){
    return Rcpp::rt(n, df) * sigma + mu;
}



// Scalar versions:

inline double d_lst(const double x, const double df, const double mu,
        const double sigma, const int log_p = 0){
    double p = R::dt((x - mu) / sigma, df, log_p);
    return log_p ? (p - log(sigma)) : ((1.0 / sigma) * p);
}

inline double p_lst(const double q, const double df, const double mu,
        const double sigma, const int lower_tail = 1, const int log_p = 0){
    return R::pt((q - mu)/sigma, df, lower_tail, log_p);
}

inline double q_lst(const double p, const double df, const double mu,
        const double sigma, const int lower_tail = 1, const int log_p = 0) {
    return R::qt(p, df, lower_tail, log_p) * sigma + mu;
}

inline double r_lst(const double df, const double mu, const double sigma){
    return R::rt(df) * sigma + mu;
}

#endif
