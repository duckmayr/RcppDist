context("Truncated Location-Scale t Distribution")

dlst <- function(x, df, mu, sigma, log_p = FALSE) {
    p = dt((x - mu)/sigma, df, log = log_p)
    if ( log_p ) {
        return(p - log(sigma))
    }
    return(p / sigma)
}

plst <- function(q, df, mu, sigma, lower_tail = TRUE, log_p = FALSE) {
    return(pt((q - mu)/sigma, df, lower.tail = lower_tail, log.p = log_p))
}

qlst <- function(p, df, mu, sigma, lower_tail = TRUE, log_p = FALSE) {
    return(qt(p, df, lower.tail = lower_tail, log.p = log_p) * sigma + mu)
}

test_that('The density functions provide correct answers', {
    dtrunclst_ <- function(x, df, mu, sigma, a, b, log_p = FALSE) {
        scale = plst(b, df, mu, sigma) - plst(a, df, mu, sigma)
        if ( log_p ) {
            result = dlst(x, df, mu, sigma, log = TRUE) - log(scale)
            return(ifelse(x < a | x > b, -Inf, result))
        }
        else {
            result = dlst(x, df, mu, sigma) / scale
            return(ifelse(x < a | x > b, 0, result))
        }
    }
    x <- seq(from = -3, to = 3)
    df <- 1; mu <- 1; sigma <- 2; a <- -2.5; b <- 2.5
    expect_equal(test_dtrunclst(x, df, mu, sigma, a, b),
                 list(
                    "VectorLog" = dtrunclst_(x, df, mu, sigma, a, b, TRUE),
                    "VectorNoLog" = dtrunclst_(x, df, mu, sigma, a, b, FALSE),
                    "DoubleLog" = dtrunclst_(x[1], df, mu, sigma, a, b, TRUE),
                    "DoubleNoLog" = dtrunclst_(x[1], df, mu, sigma, a, b, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_dtrunclst(x, df, mu, sigma, a, b),
                 list(
                    "VectorLog" = dtrunclst_(x, df, mu, sigma, a, b, TRUE),
                    "VectorNoLog" = dtrunclst_(x, df, mu, sigma, a, b, FALSE),
                    "DoubleLog" = dtrunclst_(x[1], df, mu, sigma, a, b, TRUE),
                    "DoubleNoLog" = dtrunclst_(x[1], df, mu, sigma, a, b, FALSE)
                    )
                )
})


test_that('The distribution functions provide correct answers', {
    ptrunclst_ <- function(q, df, mu, sigma, a, b, lower_tail = TRUE, log_p = FALSE) {
        scale <- plst(b, df, mu, sigma) - plst(a, df, mu, sigma)
        qq <- pmax(pmin(q, b), a)
        p <- (plst(qq, df, mu, sigma) - plst(a, df, mu, sigma)) / scale
        if ( !lower_tail ) {
            p <- 1 - p
        }
        if ( log_p ) {
            p <- log(p)
        }
        return(p)
    }
    x <- seq(from = -3, to = 3)
    df <- 1; mu <- 1; sigma <- 2; a <- -2.5; b <- 2.5
    expect_equal(test_ptrunclst(x, df, mu, sigma, a, b),
                 list(
                    "VectorLog" = ptrunclst_(x, df, mu, sigma, a, b, TRUE, TRUE),
                    "VectorNoLog" = ptrunclst_(x, df, mu, sigma, a, b),
                    "DoubleLog" = ptrunclst_(x[1], df, mu, sigma, a, b, TRUE, TRUE),
                    "DoubleNoLog" = ptrunclst_(x[1], df, mu, sigma, a, b),
                    "VectorLogNoLower" = ptrunclst_(x, df, mu, sigma, a, b, FALSE, TRUE),
                    "VectorNoLogNoLower" = ptrunclst_(x, df, mu, sigma, a, b, FALSE),
                    "DoubleLogNoLower" = ptrunclst_(x[1], df, mu, sigma, a, b, FALSE, TRUE),
                    "DoubleNoLogNoLower" = ptrunclst_(x[1], df, mu, sigma, a, b, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_ptrunclst(x, df, mu, sigma, a, b),
                 list(
                    "VectorLog" = ptrunclst_(x, df, mu, sigma, a, b, TRUE, TRUE),
                    "VectorNoLog" = ptrunclst_(x, df, mu, sigma, a, b),
                    "DoubleLog" = ptrunclst_(x[1], df, mu, sigma, a, b, TRUE, TRUE),
                    "DoubleNoLog" = ptrunclst_(x[1], df, mu, sigma, a, b),
                    "VectorLogNoLower" = ptrunclst_(x, df, mu, sigma, a, b, FALSE, TRUE),
                    "VectorNoLogNoLower" = ptrunclst_(x, df, mu, sigma, a, b, FALSE),
                    "DoubleLogNoLower" = ptrunclst_(x[1], df, mu, sigma, a, b, FALSE, TRUE),
                    "DoubleNoLogNoLower" = ptrunclst_(x[1], df, mu, sigma, a, b, FALSE)
                    )
                )
})


test_that('The quantile functions provide correct answers', {
    qtrunclst_ <- function(p, df, mu, sigma, a, b, lower_tail = TRUE, log_p = FALSE) {
        if ( log_p ) {
            p <- exp(p)
        }
        if ( !lower_tail ) {
            p <- 1 - p
        }
        F_b <- plst(b, mu, sigma, df)
        F_a <- plst(a, mu, sigma, df)
        q <- qlst(F_a + p * (F_b - F_a), df, mu, sigma)
        return(pmin(pmax(q, a), b))
    }
    x <- c(0, 0.5, 1)
    df <- 1; mu <- 1; sigma <- 2; a <- -2.5; b <- 2.5
    expect_equal(test_qtrunclst_nolog(x, df, mu, sigma, a, b),
                 list(
                    "VectorNoLog" = qtrunclst_(x, df, mu, sigma, a, b),
                    "DoubleNoLog" = qtrunclst_(x[1], df, mu, sigma, a, b),
                    "VectorNoLogNoLower" = qtrunclst_(x, df, mu, sigma, a, b, FALSE),
                    "DoubleNoLogNoLower" = qtrunclst_(x[1], df, mu, sigma, a, b, FALSE)
                    )
                )
    x <- c(-1, -2, -10)
    expect_equal(test_qtrunclst_log(x, df, mu, sigma, a, b),
                 list(
                    "VectorLog" = qtrunclst_(x, df, mu, sigma, a, b, TRUE, TRUE),
                    "DoubleLog" = qtrunclst_(x[1], df, mu, sigma, a, b, TRUE, TRUE),
                    "VectorLogNoLower" = qtrunclst_(x, df, mu, sigma, a, b, FALSE, TRUE),
                    "DoubleLogNoLower" = qtrunclst_(x[1], df, mu, sigma, a, b, FALSE, TRUE)
                    )
                )
})

