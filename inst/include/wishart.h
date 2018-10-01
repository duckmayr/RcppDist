// wishart.h
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

#ifndef RCPPDIST_WISHART_H
#define RCPPDIST_WISHART_H

inline double mvgamma(const int p, const double x) {
    double result = pow(M_PI, (p * (p - 1) * 0.25));
    for ( int j = 1; j < (p + 1); ++j ) {
        result *= R::gammafn(x - ((j - 1.0) * 0.5));
    }
    return result;
}

inline double lmvgamma(const int p, const double x) {
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

inline double diwish(const arma::mat& X, const int df, const arma::mat& S,
        const bool log_p = false) {
    double x = df * 0.5;
    int p = X.n_cols;
    if ( log_p ) {
        double P = (x * log(arma::det(S))) - (0.5 * arma::trace(S * X.i()));
        P -= ( (0.5 * (df + p + 1)) * log(arma::det(X)) );
        P -= ( (0.5 * (df * p)) * M_LN2 );
        return P - lmvgamma(p, x);
    }
    double P = pow(arma::det(S), x);
    P *= exp(-0.5 * arma::trace(S * X.i()));
    P *= pow(arma::det(X), (-1.0 * (0.5 * (df + p + 1))));
    return P / ( pow(2.0, (0.5 * (df * p))) * mvgamma(p, x) );
}

inline arma::mat riwish(const int df, const arma::mat& S) {
    return rwish(df, S.i()).i();
}

#endif
