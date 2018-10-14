#include <RcppArmadillo.h>
#include <mvnorm.h>

//' bayeslm
//'
//' Demonstrates the use of RcppDist in C++ with Bayesian linear regression
//'
//' To see an example of using RcppDist C++ functions in C++ code,
//' we can code up a Bayesian linear regression with completely uninformative
//' priors (such that estimates should be equivalent to classical estimates).
//' The code to do so is as follows:
//' \preformatted{
//' Rcpp::List bayeslm(const arma::vec& y, const arma::mat x,
//'                    const int iters = 1000) {
//'     int n = x.n_rows;
//'     int p = x.n_cols;
//'     double a = (n - p) / 2.0;
//'     arma::mat xtx = x.t() * x;
//'     arma::mat xtxinv = xtx.i();
//'     arma::vec mu = xtxinv * x.t() * y;
//'     arma::mat px = x * xtxinv * x.t();
//'     double ssq = arma::as_scalar(y.t() * (arma::eye(n, n) - px) * y);
//'     ssq *= (1.0 / (n - p));
//'     double b = 1.0 / (a * ssq);
//'     arma::mat beta_draws(iters, p);
//'     Rcpp::NumericVector sigma_draws(iters);
//'     for ( int iter = 0; iter < iters; ++iter ) {
//'         double sigmasq = 1.0 / R::rgamma(a, b);
//'         sigma_draws[iter] = sigmasq;
//'         // Here we can use our multivariate normal generator
//'         beta_draws.row(iter) = rmvnorm(1, mu, xtxinv * sigmasq);
//'     }
//'     return Rcpp::List::create(Rcpp::_["beta_draws"] = beta_draws,
//'                               Rcpp::_["sigma_draws"] = sigma_draws);
//' }
//' }
//'
//' @param y A numeric vector -- the response
//' @param x A numeric matrix -- the explanatory variables; note this assumes
//'   you have included a column of ones if you intend there to be an intercept.
//' @param iters An integer vector of length one, the number of posterior draws
//'   desired; the default is 1000.
//'
//' @return A list of length two; the first element is a numeric matrix of the
//'   beta draws and the second element is a numeric vector of the sigma draws
//' @examples
//' set.seed(123)
//' n <- 30
//' x <- cbind(1, matrix(rnorm(n*3), ncol = 3))
//' beta <- matrix(c(10, 2, -1, 3), nrow = 4)
//' y <- x %*% beta + rnorm(n)
//' freqmod <- lm(y ~ x[ , -1])
//' bayesmod <- bayeslm(y, x)
//' round(unname(coef(freqmod)), 2)
//' round(apply(bayesmod$beta_draws, 2, mean), 2)
//' c(beta)
//' @export
// [[Rcpp::export]]
Rcpp::List bayeslm(const arma::vec& y, const arma::mat x,
                   const int iters = 1000) {
    int n = x.n_rows;
    int p = x.n_cols;
    double a = (n - p) / 2.0;
    arma::mat xtx = x.t() * x;
    arma::mat xtxinv = xtx.i();
    arma::vec mu = xtxinv * x.t() * y;
    arma::mat px = x * xtxinv * x.t();
    double ssq = arma::as_scalar(y.t() * (arma::eye(n, n) - px) * y);
    ssq *= (1.0 / (n - p));
    double b = 1.0 / (a * ssq);
    arma::mat beta_draws(iters, p);
    Rcpp::NumericVector sigma_draws(iters);
    for ( int iter = 0; iter < iters; ++iter ) {
        double sigmasq = 1.0 / R::rgamma(a, b);
        sigma_draws[iter] = sigmasq;
        // Here we can use our multivariate normal generator
        beta_draws.row(iter) = rmvnorm(1, mu, xtxinv * sigmasq);
    }
    return Rcpp::List::create(Rcpp::_["beta_draws"] = beta_draws,
                              Rcpp::_["sigma_draws"] = sigma_draws);
}
