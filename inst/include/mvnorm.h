// mvnorm.h
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

#ifndef RCPPDIST_MVNORM_H
#define RCPPDIST_MVNORM_H 

inline arma::vec dmvnorm(const arma::mat& x, const arma::vec& mu,
        const arma::mat& S, const bool log_p = false) {
    arma::uword n = x.n_rows, m = x.n_cols;
    double det_S = arma::det(S);
    arma::mat S_inv = S.i();
    arma::vec result(n);
    arma::rowvec X(m);
    if ( log_p ) {
        double P = -1.0 * (x.n_cols/2.0) * M_LN_2PI - 0.5 * log(det_S);
        for ( arma::uword i = 0; i < n; ++i ) {
            X = x.row(i) - mu.t();
            result[i] = arma::as_scalar(P - 0.5 * X * S_inv * X.t());
        }
        return result;
    }
    double P = 1.0 / sqrt(pow(M_2PI, m) * det_S);
    for ( arma::uword i = 0; i < n; ++i ) {
        X = x.row(i) - mu.t();
        result[i] = arma::as_scalar(P * exp(-0.5 * X * S_inv * X.t()));
    }
    return result;
}

inline arma::mat rmvnorm(const arma::uword n, const arma::vec& mu,
        const arma::mat& S) {
    arma::uword m = S.n_cols, i, j;
    arma::mat result(n, m);
    arma::rowvec Mu = mu.t();
    for ( i = 0; i < n; ++i ) {
        for ( j = 0; j < m; ++j ) {
            result(i, j) = R::rnorm(0.0, 1.0);
        }
    }
    result = result * arma::chol(S);
    for ( i = 0; i < n; ++i ) {
        result.row(i) = result.row(i) + Mu;
    }
    return result;
}

#endif
