#include <RcppDist.h>

using namespace Rcpp;

// FOUR PARAMETER BETA DISTRIBUTION

// [[Rcpp::export]]
List test_d4beta(NumericVector x, double shape1, double shape2,
        double a, double b) {
    return List::create(
        _["VectorLog"] = d4beta(x, shape1, shape2, a, b, true),
        _["VectorNoLog"] = d4beta(x, shape1, shape2, a, b),
        _["DoubleLog"] = d_4beta(x[0], shape1, shape2, a, b, 1),
        _["DoubleNoLog"] = d_4beta(x[0], shape1, shape2, a, b)
    );
}

// [[Rcpp::export]]
List test_p4beta(NumericVector x, double shape1, double shape2,
        double a, double b) {
    return List::create(
        _["VectorLog"] = p4beta(x, shape1, shape2, a, b, true, true),
        _["VectorNoLog"] = p4beta(x, shape1, shape2, a, b),
        _["DoubleLog"] = p_4beta(x[0], shape1, shape2, a, b, 1, 1),
        _["DoubleNoLog"] = p_4beta(x[0], shape1, shape2, a, b),
        _["VectorLogNoLower"] = p4beta(x, shape1, shape2, a, b, false, true),
        _["VectorNoLogNoLower"] = p4beta(x, shape1, shape2, a, b, false),
        _["DoubleLogNoLower"] = p_4beta(x[0], shape1, shape2, a, b, 0, 1),
        _["DoubleNoLogNoLower"] = p_4beta(x[0], shape1, shape2, a, b, 0)
    );
}

// [[Rcpp::export]]
List test_q4beta_nolog(NumericVector x, double shape1, double shape2,
        double a, double b) {
    return List::create(
        _["VectorNoLog"] = q4beta(x, shape1, shape2, a, b),
        _["DoubleNoLog"] = q_4beta(x[0], shape1, shape2, a, b),
        _["VectorNoLogNoLower"] = q4beta(x, shape1, shape2, a, b, false),
        _["DoubleNoLogNoLower"] = q_4beta(x[0], shape1, shape2, a, b, 0)
    );
}

// [[Rcpp::export]]
List test_q4beta_log(NumericVector x, double shape1, double shape2,
        double a, double b) {
    return List::create(
        _["VectorLog"] = q4beta(x, shape1, shape2, a, b, true, true),
        _["DoubleLog"] = q_4beta(x[0], shape1, shape2, a, b, 1, 1),
        _["VectorLogNoLower"] = q4beta(x, shape1, shape2, a, b, false, true),
        _["DoubleLogNoLower"] = q_4beta(x[0], shape1, shape2, a, b, 0, 1)
    );
}


// LOCATION SCALE T DISTRIBUTION

// [[Rcpp::export]]
List test_dlst(NumericVector x, double df, double mu, double sigma) {
    return List::create(
        _["VectorLog"] = dlst(x, df, mu, sigma, true),
        _["VectorNoLog"] = dlst(x, df, mu, sigma),
        _["DoubleLog"] = d_lst(x[0], df, mu, sigma, 1),
        _["DoubleNoLog"] = d_lst(x[0], df, mu, sigma)
    );
}

// [[Rcpp::export]]
List test_plst(NumericVector x, double df, double mu, double sigma) {
    return List::create(
        _["VectorLog"] = plst(x, df, mu, sigma, true, true),
        _["VectorNoLog"] = plst(x, df, mu, sigma),
        _["DoubleLog"] = p_lst(x[0], df, mu, sigma, 1, 1),
        _["DoubleNoLog"] = p_lst(x[0], df, mu, sigma),
        _["VectorLogNoLower"] = plst(x, df, mu, sigma, false, true),
        _["VectorNoLogNoLower"] = plst(x, df, mu, sigma, false),
        _["DoubleLogNoLower"] = p_lst(x[0], df, mu, sigma, 0, 1),
        _["DoubleNoLogNoLower"] = p_lst(x[0], df, mu, sigma, 0)
    );
}

// [[Rcpp::export]]
List test_qlst_nolog(NumericVector x, double df, double mu, double sigma) {
    return List::create(
        _["VectorNoLog"] = qlst(x, df, mu, sigma),
        _["DoubleNoLog"] = q_lst(x[0], df, mu, sigma),
        _["VectorNoLogNoLower"] = qlst(x, df, mu, sigma, false),
        _["DoubleNoLogNoLower"] = q_lst(x[0], df, mu, sigma, 0)
    );
}

// [[Rcpp::export]]
List test_qlst_log(NumericVector x, double df, double mu, double sigma) {
    return List::create(
        _["VectorLog"] = qlst(x, df, mu, sigma, true, true),
        _["DoubleLog"] = q_lst(x[0], df, mu, sigma, 1, 1),
        _["VectorLogNoLower"] = qlst(x, df, mu, sigma, false, true),
        _["DoubleLogNoLower"] = q_lst(x[0], df, mu, sigma, 0, 1)
    );
}

// TRUNCATED NORMAL DISTRIBUTION

// [[Rcpp::export]]
List test_dtruncnorm(NumericVector x, double mu, double sigma,
        double a, double b) {
    return List::create(
        _["VectorLog"] = dtruncnorm(x, mu, sigma, a, b, true),
        _["VectorNoLog"] = dtruncnorm(x, mu, sigma, a, b),
        _["DoubleLog"] = d_truncnorm(x[0], mu, sigma, a, b, 1),
        _["DoubleNoLog"] = d_truncnorm(x[0], mu, sigma, a, b)
    );
}

// [[Rcpp::export]]
List test_ptruncnorm(NumericVector x, double mu, double sigma,
        double a, double b) {
    return List::create(
        _["VectorLog"] = ptruncnorm(x, mu, sigma, a, b, true, true),
        _["VectorNoLog"] = ptruncnorm(x, mu, sigma, a, b),
        _["DoubleLog"] = p_truncnorm(x[0], mu, sigma, a, b, 1, 1),
        _["DoubleNoLog"] = p_truncnorm(x[0], mu, sigma, a, b),
        _["VectorLogNoLower"] = ptruncnorm(x, mu, sigma, a, b, false, true),
        _["VectorNoLogNoLower"] = ptruncnorm(x, mu, sigma, a, b, false),
        _["DoubleLogNoLower"] = p_truncnorm(x[0], mu, sigma, a, b, 0, 1),
        _["DoubleNoLogNoLower"] = p_truncnorm(x[0], mu, sigma, a, b, 0)
    );
}

// [[Rcpp::export]]
List test_qtruncnorm_nolog(NumericVector x, double mu, double sigma,
        double a, double b) {
    return List::create(
        _["VectorNoLog"] = qtruncnorm(x, mu, sigma, a, b),
        _["DoubleNoLog"] = q_truncnorm(x[0], mu, sigma, a, b),
        _["VectorNoLogNoLower"] = qtruncnorm(x, mu, sigma, a, b, false),
        _["DoubleNoLogNoLower"] = q_truncnorm(x[0], mu, sigma, a, b, 0)
    );
}

// [[Rcpp::export]]
List test_qtruncnorm_log(NumericVector x, double mu, double sigma,
        double a, double b) {
    return List::create(
        _["VectorLog"] = qtruncnorm(x, mu, sigma, a, b, true, true),
        _["DoubleLog"] = q_truncnorm(x[0], mu, sigma, a, b, 1, 1),
        _["VectorLogNoLower"] = qtruncnorm(x, mu, sigma, a, b, false, true),
        _["DoubleLogNoLower"] = q_truncnorm(x[0], mu, sigma, a, b, 0, 1)
    );
}
