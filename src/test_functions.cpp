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
} // nocov

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
} // nocov

// [[Rcpp::export]]
List test_q4beta_nolog(NumericVector x, double shape1, double shape2,
        double a, double b) {
    return List::create(
        _["VectorNoLog"] = q4beta(x, shape1, shape2, a, b),
        _["DoubleNoLog"] = q_4beta(x[0], shape1, shape2, a, b),
        _["VectorNoLogNoLower"] = q4beta(x, shape1, shape2, a, b, false),
        _["DoubleNoLogNoLower"] = q_4beta(x[0], shape1, shape2, a, b, 0)
    );
} // nocov

// [[Rcpp::export]]
List test_q4beta_log(NumericVector x, double shape1, double shape2,
        double a, double b) {
    return List::create(
        _["VectorLog"] = q4beta(x, shape1, shape2, a, b, true, true),
        _["DoubleLog"] = q_4beta(x[0], shape1, shape2, a, b, 1, 1),
        _["VectorLogNoLower"] = q4beta(x, shape1, shape2, a, b, false, true),
        _["DoubleLogNoLower"] = q_4beta(x[0], shape1, shape2, a, b, 0, 1)
    );
} // nocov



// LOCATION-SCALE T DISTRIBUTION

// [[Rcpp::export]]
List test_dlst(NumericVector x, double df, double mu, double sigma) {
    return List::create(
        _["VectorLog"] = dlst(x, df, mu, sigma, true),
        _["VectorNoLog"] = dlst(x, df, mu, sigma),
        _["DoubleLog"] = d_lst(x[0], df, mu, sigma, 1),
        _["DoubleNoLog"] = d_lst(x[0], df, mu, sigma)
    );
} // nocov

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
} // nocov

// [[Rcpp::export]]
List test_qlst_nolog(NumericVector x, double df, double mu, double sigma) {
    return List::create(
        _["VectorNoLog"] = qlst(x, df, mu, sigma),
        _["DoubleNoLog"] = q_lst(x[0], df, mu, sigma),
        _["VectorNoLogNoLower"] = qlst(x, df, mu, sigma, false),
        _["DoubleNoLogNoLower"] = q_lst(x[0], df, mu, sigma, 0)
    );
} // nocov

// [[Rcpp::export]]
List test_qlst_log(NumericVector x, double df, double mu, double sigma) {
    return List::create(
        _["VectorLog"] = qlst(x, df, mu, sigma, true, true),
        _["DoubleLog"] = q_lst(x[0], df, mu, sigma, 1, 1),
        _["VectorLogNoLower"] = qlst(x, df, mu, sigma, false, true),
        _["DoubleLogNoLower"] = q_lst(x[0], df, mu, sigma, 0, 1)
    );
} // nocov



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
} // nocov

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
} // nocov

// [[Rcpp::export]]
List test_qtruncnorm_nolog(NumericVector x, double mu, double sigma,
        double a, double b) {
    return List::create(
        _["VectorNoLog"] = qtruncnorm(x, mu, sigma, a, b),
        _["DoubleNoLog"] = q_truncnorm(x[0], mu, sigma, a, b),
        _["VectorNoLogNoLower"] = qtruncnorm(x, mu, sigma, a, b, false),
        _["DoubleNoLogNoLower"] = q_truncnorm(x[0], mu, sigma, a, b, 0)
    );
} // nocov

// [[Rcpp::export]]
List test_qtruncnorm_log(NumericVector x, double mu, double sigma,
        double a, double b) {
    return List::create(
        _["VectorLog"] = qtruncnorm(x, mu, sigma, a, b, true, true),
        _["DoubleLog"] = q_truncnorm(x[0], mu, sigma, a, b, 1, 1),
        _["VectorLogNoLower"] = qtruncnorm(x, mu, sigma, a, b, false, true),
        _["DoubleLogNoLower"] = q_truncnorm(x[0], mu, sigma, a, b, 0, 1)
    );
} // nocov



// TRUNCATED T DISTRIBUTION

// [[Rcpp::export]]
List test_dtrunct(NumericVector x, double df, double a, double b) {
    return List::create(
        _["VectorLog"] = dtrunct(x, df, a, b, true),
        _["VectorNoLog"] = dtrunct(x, df, a, b),
        _["DoubleLog"] = d_trunct(x[0], df, a, b, 1),
        _["DoubleNoLog"] = d_trunct(x[0], df, a, b)
    );
} // nocov

// [[Rcpp::export]]
List test_ptrunct(NumericVector x, double df, double a, double b) {
    return List::create(
        _["VectorLog"] = ptrunct(x, df, a, b, true, true),
        _["VectorNoLog"] = ptrunct(x, df, a, b),
        _["DoubleLog"] = p_trunct(x[0], df, a, b, 1, 1),
        _["DoubleNoLog"] = p_trunct(x[0], df, a, b),
        _["VectorLogNoLower"] = ptrunct(x, df, a, b, false, true),
        _["VectorNoLogNoLower"] = ptrunct(x, df, a, b, false),
        _["DoubleLogNoLower"] = p_trunct(x[0], df, a, b, 0, 1),
        _["DoubleNoLogNoLower"] = p_trunct(x[0], df, a, b, 0)
    );
} // nocov

// [[Rcpp::export]]
List test_qtrunct_nolog(NumericVector x, double df, double a, double b) {
    return List::create(
        _["VectorNoLog"] = qtrunct(x, df, a, b),
        _["DoubleNoLog"] = q_trunct(x[0], df, a, b),
        _["VectorNoLogNoLower"] = qtrunct(x, df, a, b, false),
        _["DoubleNoLogNoLower"] = q_trunct(x[0], df, a, b, 0)
    );
} // nocov

// [[Rcpp::export]]
List test_qtrunct_log(NumericVector x, double df, double a, double b) {
    return List::create(
        _["VectorLog"] = qtrunct(x, df, a, b, true, true),
        _["DoubleLog"] = q_trunct(x[0], df, a, b, 1, 1),
        _["VectorLogNoLower"] = qtrunct(x, df, a, b, false, true),
        _["DoubleLogNoLower"] = q_trunct(x[0], df, a, b, 0, 1)
    );
} // nocov



// TRUNCATED LOCATION-SCALE T DISTRIBUTION

// [[Rcpp::export]]
List test_dtrunclst(NumericVector x, double df, double mu, double sigma,
        double a, double b) {
    return List::create(
        _["VectorLog"] = dtrunclst(x, df, mu, sigma, a, b, true),
        _["VectorNoLog"] = dtrunclst(x, df, mu, sigma, a, b),
        _["DoubleLog"] = d_trunclst(x[0], df, mu, sigma, a, b, 1),
        _["DoubleNoLog"] = d_trunclst(x[0], df, mu, sigma, a, b)
    );
} // nocov

// [[Rcpp::export]]
List test_ptrunclst(NumericVector x, double df, double mu, double sigma,
        double a, double b) {
    return List::create(
        _["VectorLog"] = ptrunclst(x, df, mu, sigma, a, b, true, true),
        _["VectorNoLog"] = ptrunclst(x, df, mu, sigma, a, b),
        _["DoubleLog"] = p_trunclst(x[0], df, mu, sigma, a, b, 1, 1),
        _["DoubleNoLog"] = p_trunclst(x[0], df, mu, sigma, a, b),
        _["VectorLogNoLower"] = ptrunclst(x, df, mu, sigma, a, b, false, true),
        _["VectorNoLogNoLower"] = ptrunclst(x, df, mu, sigma, a, b, false),
        _["DoubleLogNoLower"] = p_trunclst(x[0], df, mu, sigma, a, b, 0, 1),
        _["DoubleNoLogNoLower"] = p_trunclst(x[0], df, mu, sigma, a, b, 0)
    );
} // nocov

// [[Rcpp::export]]
List test_qtrunclst_nolog(NumericVector x, double df, double mu, double sigma,
        double a, double b) {
    return List::create(
        _["VectorNoLog"] = qtrunclst(x, df, mu, sigma, a, b),
        _["DoubleNoLog"] = q_trunclst(x[0], df, mu, sigma, a, b),
        _["VectorNoLogNoLower"] = qtrunclst(x, df, mu, sigma, a, b, false),
        _["DoubleNoLogNoLower"] = q_trunclst(x[0], df, mu, sigma, a, b, 0)
    );
} // nocov

// [[Rcpp::export]]
List test_qtrunclst_log(NumericVector x, double df, double mu, double sigma,
        double a, double b) {
    return List::create(
        _["VectorLog"] = qtrunclst(x, df, mu, sigma, a, b, true, true),
        _["DoubleLog"] = q_trunclst(x[0], df, mu, sigma, a, b, 1, 1),
        _["VectorLogNoLower"] = qtrunclst(x, df, mu, sigma, a, b, false, true),
        _["DoubleLogNoLower"] = q_trunclst(x[0], df, mu, sigma, a, b, 0, 1)
    );
} // nocov



// FOUR PARAMETER BETA DISTRIBUTION

// [[Rcpp::export]]
List test_dtri(NumericVector x, double a, double b, double c) {
    return List::create(
        _["VectorLog"] = dtri(x, a, b, c, true),
        _["VectorNoLog"] = dtri(x, a, b, c),
        _["DoubleLog"] = d_tri(x[0], a, b, c, 1),
        _["DoubleNoLog"] = d_tri(x[0], a, b, c)
    );
} // nocov

// [[Rcpp::export]]
List test_ptri(NumericVector x, double a, double b, double c) {
    return List::create(
        _["VectorLog"] = ptri(x, a, b, c, true, true),
        _["VectorNoLog"] = ptri(x, a, b, c),
        _["DoubleLog"] = p_tri(x[0], a, b, c, 1, 1),
        _["DoubleNoLog"] = p_tri(x[0], a, b, c),
        _["VectorLogNoLower"] = ptri(x, a, b, c, false, true),
        _["VectorNoLogNoLower"] = ptri(x, a, b, c, false),
        _["DoubleLogNoLower"] = p_tri(x[0], a, b, c, 0, 1),
        _["DoubleNoLogNoLower"] = p_tri(x[0], a, b, c, 0)
    );
} // nocov

// [[Rcpp::export]]
List test_qtri_nolog(NumericVector x, double a, double b, double c) {
    return List::create(
        _["VectorNoLog"] = qtri(x, a, b, c),
        _["DoubleNoLog"] = q_tri(x[0], a, b, c),
        _["VectorNoLogNoLower"] = qtri(x, a, b, c, false),
        _["DoubleNoLogNoLower"] = q_tri(x[0], a, b, c, 0)
    );
} // nocov

// [[Rcpp::export]]
List test_qtri_log(NumericVector x, double a, double b, double c) {
    return List::create(
        _["VectorLog"] = qtri(x, a, b, c, true, true),
        _["DoubleLog"] = q_tri(x[0], a, b, c, 1, 1),
        _["VectorLogNoLower"] = qtri(x, a, b, c, false, true),
        _["DoubleLogNoLower"] = q_tri(x[0], a, b, c, 0, 1)
    );
} // nocov



// MULTIVARIATE NORMAL DISTRIBUTION

// [[Rcpp::export]]
List test_dmvnorm(arma::mat x, arma::vec mu, arma::mat S) {
    return List::create(
        _["Log"] = dmvnorm(x, mu, S, true),
        _["NoLog"] = dmvnorm(x, mu, S, false)
    );
} // nocov



// MULTIVARIATE T DISTRIBUTION

// [[Rcpp::export]]
List test_dmvt(arma::mat x, arma::vec mu, arma::mat S, double df) {
    return List::create(
        _["Log"] = dmvt(x, mu, S, df, true),
        _["NoLog"] = dmvt(x, mu, S, df, false)
    );
} // nocov



// WISHART DISTRIBUTION

// [[Rcpp::export]]
List test_dwish(arma::mat x, int df, arma::mat S) {
    return List::create(
        _["Log"] = dwish(x, df, S, true),
        _["NoLog"] = dwish(x, df, S, false)
    );
} // nocov



// INVERSE WISHART DISTRIBUTION

// [[Rcpp::export]]
List test_diwish(arma::mat x, int df, arma::mat S) {
    return List::create(
        _["Log"] = diwish(x, df, S, true),
        _["NoLog"] = diwish(x, df, S, false)
    );
} // nocov
