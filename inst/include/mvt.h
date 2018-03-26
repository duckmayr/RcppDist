#ifndef RCPPDIST_MVT_H
#define RCPPDIST_MVT_H

#include <mvnorm.h>

inline arma::vec dmvt(const arma::mat& x, const arma::vec& mu, arma::mat S,
        const double df, const bool log_p = false) {
    arma::uword n = x.n_rows, m = x.n_cols, i;
    double det_S = arma::det(S);
    S = S.i();
    arma::vec result(n);
    arma::rowvec X(m);
    if ( log_p ) {
        double P = R::lgammafn((df + m) * 0.5) - R::lgammafn(df * 0.5);
        P -= ( (m * 0.5) * (log(df) + log(M_PI)) + 0.5 * log(det_S) );
        for ( arma::uword i = 0; i < n; ++i ) {
            X = x.row(i) - mu.t();
            result[i] = arma::as_scalar(P - ((df + m) * 0.5) * log(1.0 + (1.0 / df) * X * S * X.t()));
        }
        return result;
    }
    double P = R::gammafn((df + m) * 0.5);
    P /= (R::gammafn(df*0.5) * pow(df, m*0.5) * pow(M_PI, m*0.5) * sqrt(det_S));
    for ( arma::uword i = 0; i < n; ++i ) {
        X = x.row(i) - mu.t();
        result[i] = arma::as_scalar(P/pow(1.0+(1.0/df)*X*S*X.t(), (df+m)*0.5));
    }
    return result;
}

inline arma::mat rmvt(const int n, const arma::vec& mu, const arma::mat& S,
               const double df) {
    arma::uword m = S.n_cols;
    arma::vec U = Rcpp::rchisq(n, df);
    U = sqrt(df / U);
    arma::mat Y = rmvnorm(n, arma::vec(m, arma::fill::zeros), S).t();
    arma::mat result(m, n);
    for ( arma::uword i = 0; i < n; ++i ) {
        result.col(i) = mu + Y.col(i) * U[i];
    }
    return result.t();
}

#endif

