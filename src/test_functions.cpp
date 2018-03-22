#include <RcppDist.h>

using namespace Rcpp;

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
