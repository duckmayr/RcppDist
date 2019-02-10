context("bayeslm")

test_that('Example estimates are as expected', {
    set.seed(123)
    n <- 30
    x <- cbind(1, matrix(rnorm(n*3), ncol = 3))
    beta <- matrix(c(10, -2, -1, 3), nrow = 4)
    sigma <- 1.5
    y <- x %*% beta + rnorm(n, mean = 0, sd = sigma)
    freqmod <- lm(y ~ x[ , -1])
    bayesmod <- bayeslm(y, x, 10000)
    beta_means <- apply(bayesmod$beta_draws, 2, mean)
    sigmasq_mean <- mean(bayesmod$sigmasq_draws)
    expect_equivalent(beta_means, coef(freqmod), tolerance = 0.01)
    expect_equivalent(sqrt(sigmasq_mean), sigma(freqmod), tolerance = 0.06)
})
