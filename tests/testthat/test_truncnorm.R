context("Truncated Normal Distribution")

test_that('The density functions provide correct answers', {
    dtruncnorm_ <- function(x, mu, sigma, a, b, log_p = FALSE) {
        scale = pnorm(b, mu, sigma) - pnorm(a, mu, sigma)
        if ( log_p ) {
            result = dnorm(x, mu, sigma, log = TRUE) - log(scale)
            return(ifelse(x < a | x > b, -Inf, result))
        }
        else {
            result = dnorm(x, mu, sigma) / scale
            return(ifelse(x < a | x > b, 0, result))
        }
    }
    x <- seq(from = -3, to = 3)
    mu <- 1; sigma <- 2; a <- -2.5; b <- 2.5
    expect_equal(test_dtruncnorm(x, mu, sigma, a, b),
                 list(
                    "VectorLog" = dtruncnorm_(x, mu, sigma, a, b, TRUE),
                    "VectorNoLog" = dtruncnorm_(x, mu, sigma, a, b, FALSE),
                    "DoubleLog" = dtruncnorm_(x[1], mu, sigma, a, b, TRUE),
                    "DoubleNoLog" = dtruncnorm_(x[1], mu, sigma, a, b, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_dtruncnorm(x, mu, sigma, a, b),
                 list(
                    "VectorLog" = dtruncnorm_(x, mu, sigma, a, b, TRUE),
                    "VectorNoLog" = dtruncnorm_(x, mu, sigma, a, b, FALSE),
                    "DoubleLog" = dtruncnorm_(x[1], mu, sigma, a, b, TRUE),
                    "DoubleNoLog" = dtruncnorm_(x[1], mu, sigma, a, b, FALSE)
                    )
                )
})


test_that('The distribution functions provide correct answers', {
    ptruncnorm_ <- function(q, mu, sigma, a, b, lower_tail = TRUE,
                            log_p = FALSE) {
        scale <- pnorm(b, mu, sigma) - pnorm(a, mu, sigma)
        qq <- pmax(pmin(q, b), a)
        p <- (pnorm(qq, mu, sigma) - pnorm(a, mu, sigma)) / scale
        if ( !lower_tail ) {
            p <- 1 - p
        }
        if ( log_p ) {
            p <- log(p)
        }
        return(p)
    }
    x <- seq(from = -3, to = 3)
    mu <- 1; sigma <- 2; a <- -2.5; b <- 2.5
    expect_equal(test_ptruncnorm(x, mu, sigma, a, b),
                 list(
                    "VectorLog" = ptruncnorm_(x, mu, sigma, a, b, TRUE, TRUE),
                    "VectorNoLog" = ptruncnorm_(x, mu, sigma, a, b),
                    "DoubleLog" = ptruncnorm_(x[1], mu, sigma, a, b, TRUE, TRUE),
                    "DoubleNoLog" = ptruncnorm_(x[1], mu, sigma, a, b),
                    "VectorLogNoLower" = ptruncnorm_(x, mu, sigma, a, b, FALSE, TRUE),
                    "VectorNoLogNoLower" = ptruncnorm_(x, mu, sigma, a, b, FALSE),
                    "DoubleLogNoLower" = ptruncnorm_(x[1], mu, sigma, a, b, FALSE, TRUE),
                    "DoubleNoLogNoLower" = ptruncnorm_(x[1], mu, sigma, a, b, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_ptruncnorm(x, mu, sigma, a, b),
                 list(
                    "VectorLog" = ptruncnorm_(x, mu, sigma, a, b, TRUE, TRUE),
                    "VectorNoLog" = ptruncnorm_(x, mu, sigma, a, b),
                    "DoubleLog" = ptruncnorm_(x[1], mu, sigma, a, b, TRUE, TRUE),
                    "DoubleNoLog" = ptruncnorm_(x[1], mu, sigma, a, b),
                    "VectorLogNoLower" = ptruncnorm_(x, mu, sigma, a, b, FALSE, TRUE),
                    "VectorNoLogNoLower" = ptruncnorm_(x, mu, sigma, a, b, FALSE),
                    "DoubleLogNoLower" = ptruncnorm_(x[1], mu, sigma, a, b, FALSE, TRUE),
                    "DoubleNoLogNoLower" = ptruncnorm_(x[1], mu, sigma, a, b, FALSE)
                    )
                )
})


test_that('The quantile functions provide correct answers', {
    qtruncnorm_ <- function(p, mu, sigma, a, b, lower_tail = TRUE,
                            log_p = FALSE) {
        if ( log_p ) {
            p <- exp(p)
        }
        if ( !lower_tail ) {
            p <- 1 - p
        }
        F_b <- pnorm(b, mu, sigma)
        F_a <- pnorm(a, mu, sigma)
        q <- qnorm(F_a + p * (F_b - F_a), mu, sigma)
        return(pmin(pmax(q, a), b))
    }
    x <- c(0, 0.5, 1)
    mu <- 1; sigma <- 2; a <- -2.5; b <- 2.5
    expect_equal(test_qtruncnorm_nolog(x, mu, sigma, a, b),
                 list(
                    "VectorNoLog" = qtruncnorm_(x, mu, sigma, a, b),
                    "DoubleNoLog" = qtruncnorm_(x[1], mu, sigma, a, b),
                    "VectorNoLogNoLower" = qtruncnorm_(x, mu, sigma, a, b, FALSE),
                    "DoubleNoLogNoLower" = qtruncnorm_(x[1], mu, sigma, a, b, FALSE)
                    )
                )
    x <- c(-1, -2, -10)
    expect_equal(test_qtruncnorm_log(x, mu, sigma, a, b),
                 list(
                    "VectorLog" = qtruncnorm_(x, mu, sigma, a, b, TRUE, TRUE),
                    "DoubleLog" = qtruncnorm_(x[1], mu, sigma, a, b, TRUE, TRUE),
                    "VectorLogNoLower" = qtruncnorm_(x, mu, sigma, a, b, FALSE, TRUE),
                    "DoubleLogNoLower" = qtruncnorm_(x[1], mu, sigma, a, b, FALSE, TRUE)
                    )
                )
})

