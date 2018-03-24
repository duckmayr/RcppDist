context("Triangular Distribution")

test_that('The density functions provide correct answers', {
    dtri_ <- function(x, a, b, c, log_p = FALSE) {
        p <- ifelse(x < c, (2 * (x - a)) / ((b - a) * (c - a)),
                    (2 * (b - x)) / ((b - a) * (b - c)))
        p[x < a | x > b] <- 0
        if ( log_p ) {
            return(log(p))
        }
        return(p)
    }
    x <- seq(from = -3, to = 3)
    a <- -2.5; b <- 2.5; c <- 0;
    expect_equal(test_dtri(x, a, b, c),
                 list(
                    "VectorLog" = dtri_(x, a, b, c, TRUE),
                    "VectorNoLog" = dtri_(x, a, b, c, FALSE),
                    "DoubleLog" = dtri_(x[1], a, b, c, TRUE),
                    "DoubleNoLog" = dtri_(x[1], a, b, c, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_dtri(x, a, b, c),
                 list(
                    "VectorLog" = dtri_(x, a, b, c, TRUE),
                    "VectorNoLog" = dtri_(x, a, b, c, FALSE),
                    "DoubleLog" = dtri_(x[1], a, b, c, TRUE),
                    "DoubleNoLog" = dtri_(x[1], a, b, c, FALSE)
                    )
                )
})


test_that('The distribution functions provide correct answers', {
    ptri_ <- function(q, a, b, c, lower_tail = TRUE, log_p = FALSE) {
        p <- ifelse(q < c, ((x - a)^2) / ((b - a) * (c - a)),
                    1 - ((x - b)^2) / ((b - a) * (b - c)))
        p[q < a] <- 0
        p[q > b] <- 1
        if ( log_p ) {
            if ( lower_tail ) {
                return(log(p))
            }
            return(log(1 - p))
        }
        if ( lower_tail ) {
            return(p)
        }
        return(1 - p)
    }
    x <- seq(from = -3, to = 3)
    a <- -2.5; b <- 2.5; c <- 0;
    expect_equal(test_ptri(x, a, b, c),
                 list(
                    "VectorLog" = ptri_(x, a, b, c, TRUE, TRUE),
                    "VectorNoLog" = ptri_(x, a, b, c),
                    "DoubleLog" = ptri_(x[1], a, b, c, TRUE, TRUE),
                    "DoubleNoLog" = ptri_(x[1], a, b, c),
                    "VectorLogNoLower" = ptri_(x, a, b, c, FALSE, TRUE),
                    "VectorNoLogNoLower" = ptri_(x, a, b, c, FALSE),
                    "DoubleLogNoLower" = ptri_(x[1], a, b, c, FALSE, TRUE),
                    "DoubleNoLogNoLower" = ptri_(x[1], a, b, c, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_ptri(x, a, b, c),
                 list(
                    "VectorLog" = ptri_(x, a, b, c, TRUE, TRUE),
                    "VectorNoLog" = ptri_(x, a, b, c),
                    "DoubleLog" = ptri_(x[1], a, b, c, TRUE, TRUE),
                    "DoubleNoLog" = ptri_(x[1], a, b, c),
                    "VectorLogNoLower" = ptri_(x, a, b, c, FALSE, TRUE),
                    "VectorNoLogNoLower" = ptri_(x, a, b, c, FALSE),
                    "DoubleLogNoLower" = ptri_(x[1], a, b, c, FALSE, TRUE),
                    "DoubleNoLogNoLower" = ptri_(x[1], a, b, c, FALSE)
                    )
                )
})


test_that('The quantile functions provide correct answers', {
    qtri_ <- function(p, a, b, c, lower_tail = TRUE, log_p = FALSE) {
        if ( log_p ) {
            p <- exp(p)
        }
        if ( !lower_tail ) {
            p <- 1 - p
        }
        cutoff <- (c - a) / (b - a)
        q <- numeric(length(p))
        q[p < cutoff] <- a + sqrt((b - a) * (c - a) * p[p < cutoff])
        q[p >= cutoff] <- b - sqrt((b - a) * (b - c) * (1 - p[p >= cutoff]))
        return(q)
    }
    x <- c(0, 0.5, 1)
    a <- -2.5; b <- 2.5; c <- 0;
    expect_equal(test_qtri_nolog(x, a, b, c),
                 list(
                    "VectorNoLog" = qtri_(x, a, b, c),
                    "DoubleNoLog" = qtri_(x[1], a, b, c),
                    "VectorNoLogNoLower" = qtri_(x, a, b, c, FALSE),
                    "DoubleNoLogNoLower" = qtri_(x[1], a, b, c, FALSE)
                    )
                )
    x <- c(0, -1, -2, -10)
    expect_equal(test_qtri_log(x, a, b, c),
                 list(
                    "VectorLog" = qtri_(x, a, b, c, TRUE, TRUE),
                    "DoubleLog" = qtri_(x[1], a, b, c, TRUE, TRUE),
                    "VectorLogNoLower" = qtri_(x, a, b, c, FALSE, TRUE),
                    "DoubleLogNoLower" = qtri_(x[1], a, b, c, FALSE, TRUE)
                    )
                )
})

