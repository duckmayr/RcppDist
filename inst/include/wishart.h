#ifndef RCPPDIST_WISHART_H
#define RCPPDIST_WISHART_H

inline double mvgamma(const int p, const int x) {
    double result = pow(M_PI, (p * (p - 1) * 0.25));
    for ( int j = 1; j < (p + 1); ++j ) {
        result *= R::gammafn(x - ((j - 1.0) * 0.5));
    }
    return result;
}

inline double lmvgamma(const int p, const int x) {
    double result = (p * (p - 1) * 0.25) * log(M_PI);
    for ( int j = 1; j < (p + 1); ++j ) {
        result += R::lgammafn(x - ((j - 1.0) * 0.5));
    }
    return result;
}

inline double dwish(const arma::mat& X, const int df, const arma::mat& S,
        const bool log_p = false) {
    double x = df * 0.5;
    int p = X.n_cols;
    if ( log_p ) {
        double P = ((df - p - 1) * 0.5) * log(arma::det(X));
        P -= (arma::trace(S.i() * X) * 0.5);
        P -= ( (p * x * M_LN2) + (x * log(arma::det(S))) );
        return P - lmvgamma(p, x);
    }
    double P = pow(arma::det(X), ((df - p - 1.0) * 0.5));
    P *= exp(-0.5 * arma::trace(S.i() * X));
    return P / ( pow(2.0, (x * p)) * pow(arma::det(S), x) * mvgamma(p, x) );
}

inline arma::mat rwish(const int df, const arma::mat& S) {
    arma::uword m = S.n_cols;
    arma::uword i, j;
    arma::mat A(m, m, arma::fill::zeros);
    for ( i = 1; i < m; ++i ) {
        for ( j = 0; j < i; ++j ) {
            A(i, j) = R::rnorm(0.0, 1.0);
        }
    }
    for ( i = 0; i < m; ++i ) {
        A(i, i) = sqrt(R::rchisq(df - i));
    }
    arma::mat B = A.t() * arma::chol(S);
    return B.t() * B;
}

#endif
