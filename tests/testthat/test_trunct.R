context("Truncated t Distribution")

test_that('The density functions provide correct answers', {
    dtrunct_ <- function(x, df, a, b, log_p = FALSE) {
        scale = pt(b, df) - pt(a, df)
        if ( log_p ) {
            result = dt(x, df, log = TRUE) - log(scale)
            return(ifelse(x < a | x > b, -Inf, result))
        }
        else {
            result = dt(x, df) / scale
            return(ifelse(x < a | x > b, 0, result))
        }
    }
    x <- seq(from = -3, to = 3)
    df <- 1; a <- -2.5; b <- 2.5
    expect_equal(test_dtrunct(x, df, a, b),
                 list(
                    "VectorLog" = dtrunct_(x, df, a, b, TRUE),
                    "VectorNoLog" = dtrunct_(x, df, a, b, FALSE),
                    "DoubleLog" = dtrunct_(x[1], df, a, b, TRUE),
                    "DoubleNoLog" = dtrunct_(x[1], df, a, b, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_dtrunct(x, df, a, b),
                 list(
                    "VectorLog" = dtrunct_(x, df, a, b, TRUE),
                    "VectorNoLog" = dtrunct_(x, df, a, b, FALSE),
                    "DoubleLog" = dtrunct_(x[1], df, a, b, TRUE),
                    "DoubleNoLog" = dtrunct_(x[1], df, a, b, FALSE)
                    )
                )
})


test_that('The distribution functions provide correct answers', {
    ptrunct_ <- function(q, df, a, b, lower_tail = TRUE, log_p = FALSE) {
        scale <- pt(b, df) - pt(a, df)
        qq <- pmax(pmin(q, b), a)
        p <- (pt(qq, df) - pt(a, df)) / scale
        if ( !lower_tail ) {
            p <- 1 - p
        }
        if ( log_p ) {
            p <- log(p)
        }
        return(p)
    }
    x <- seq(from = -3, to = 3)
    df <- 1; a <- -2.5; b <- 2.5
    expect_equal(test_ptrunct(x, df, a, b),
                 list(
                    "VectorLog" = ptrunct_(x, df, a, b, TRUE, TRUE),
                    "VectorNoLog" = ptrunct_(x, df, a, b),
                    "DoubleLog" = ptrunct_(x[1], df, a, b, TRUE, TRUE),
                    "DoubleNoLog" = ptrunct_(x[1], df, a, b),
                    "VectorLogNoLower" = ptrunct_(x, df, a, b, FALSE, TRUE),
                    "VectorNoLogNoLower" = ptrunct_(x, df, a, b, FALSE),
                    "DoubleLogNoLower" = ptrunct_(x[1], df, a, b, FALSE, TRUE),
                    "DoubleNoLogNoLower" = ptrunct_(x[1], df, a, b, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_ptrunct(x, df, a, b),
                 list(
                    "VectorLog" = ptrunct_(x, df, a, b, TRUE, TRUE),
                    "VectorNoLog" = ptrunct_(x, df, a, b),
                    "DoubleLog" = ptrunct_(x[1], df, a, b, TRUE, TRUE),
                    "DoubleNoLog" = ptrunct_(x[1], df, a, b),
                    "VectorLogNoLower" = ptrunct_(x, df, a, b, FALSE, TRUE),
                    "VectorNoLogNoLower" = ptrunct_(x, df, a, b, FALSE),
                    "DoubleLogNoLower" = ptrunct_(x[1], df, a, b, FALSE, TRUE),
                    "DoubleNoLogNoLower" = ptrunct_(x[1], df, a, b, FALSE)
                    )
                )
})


test_that('The quantile functions provide correct answers', {
    qtrunct_ <- function(p, df, a, b, lower_tail = TRUE, log_p = FALSE) {
        if ( log_p ) {
            p <- exp(p)
        }
        if ( !lower_tail ) {
            p <- 1 - p
        }
        F_b <- pt(b, df)
        F_a <- pt(a, df)
        q <- qt(F_a + p * (F_b - F_a), df)
        return(pmin(pmax(q, a), b))
    }
    x <- c(0, 0.5, 1)
    df <- 1; a <- -2.5; b <- 2.5
    expect_equal(test_qtrunct_nolog(x, df, a, b),
                 list(
                    "VectorNoLog" = qtrunct_(x, df, a, b),
                    "DoubleNoLog" = qtrunct_(x[1], df, a, b),
                    "VectorNoLogNoLower" = qtrunct_(x, df, a, b, FALSE),
                    "DoubleNoLogNoLower" = qtrunct_(x[1], df, a, b, FALSE)
                    )
                )
    x <- c(-1, -2, -10)
    expect_equal(test_qtrunct_log(x, df, a, b),
                 list(
                    "VectorLog" = qtrunct_(x, df, a, b, TRUE, TRUE),
                    "DoubleLog" = qtrunct_(x[1], df, a, b, TRUE, TRUE),
                    "VectorLogNoLower" = qtrunct_(x, df, a, b, FALSE, TRUE),
                    "DoubleLogNoLower" = qtrunct_(x[1], df, a, b, FALSE, TRUE)
                    )
                )
})

