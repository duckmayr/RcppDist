context("Four Parameter Beta Distribution")

test_that('The density functions provide correct answers', {
    d4beta_ <- function(x, shape1, shape2, a, b, log_p = FALSE) {
        res <- dbeta((x-a) / (b-a), shape1, shape2, log = log_p)
        if ( log_p ) {
            res <- res - log(b - a)
            res <- ifelse(x < a | x > b, -Inf, res)
        }
        else {
            res <- res / (b - a)
            res <- ifelse(x < a | x > b, 0, res)
        }
        return(res)
    }
    x <- c(-3, 2, 0, 4, -1)
    s1 <- 2.0; s2 <- 2.0; a <- -2.5; b <- 2.5;
    expect_equal(test_d4beta(x, s1, s2, a, b),
                 list(
                    "VectorLog" = d4beta_(x, s1, s2, a, b, TRUE),
                    "VectorNoLog" = d4beta_(x, s1, s2, a, b, FALSE),
                    "DoubleLog" = d4beta_(x[1], s1, s2, a, b, TRUE),
                    "DoubleNoLog" = d4beta_(x[1], s1, s2, a, b, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_d4beta(x, s1, s2, a, b),
                 list(
                    "VectorLog" = d4beta_(x, s1, s2, a, b, TRUE),
                    "VectorNoLog" = d4beta_(x, s1, s2, a, b, FALSE),
                    "DoubleLog" = d4beta_(x[1], s1, s2, a, b, TRUE),
                    "DoubleNoLog" = d4beta_(x[1], s1, s2, a, b, FALSE)
                    )
                )
})


test_that('The distribution functions provide correct answers', {
    p4beta_ <- function(q, shape1, shape2, a, b,
                        lower_tail = TRUE, log_p = FALSE) {
        return(pbeta((q - a) / (b - a), shape1, shape2,
                     lower.tail = lower_tail, log.p = log_p))
    }
    x <- c(-3, 2, 0, 4, -1)
    s1 <- 2.0; s2 <- 2.0; a <- -2.5; b <- 2.5;
    expect_equal(test_p4beta(x, s1, s2, a, b),
                 list(
                    "VectorLog" = p4beta_(x, s1, s2, a, b, TRUE, TRUE),
                    "VectorNoLog" = p4beta_(x, s1, s2, a, b),
                    "DoubleLog" = p4beta_(x[1], s1, s2, a, b, TRUE, TRUE),
                    "DoubleNoLog" = p4beta_(x[1], s1, s2, a, b),
                    "VectorLogNoLower" = p4beta_(x, s1, s2, a, b, FALSE, TRUE),
                    "VectorNoLogNoLower" = p4beta_(x, s1, s2, a, b, FALSE),
                    "DoubleLogNoLower" = p4beta_(x[1], s1, s2, a, b, FALSE, TRUE),
                    "DoubleNoLogNoLower" = p4beta_(x[1], s1, s2, a, b, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_p4beta(x, s1, s2, a, b),
                 list(
                    "VectorLog" = p4beta_(x, s1, s2, a, b, TRUE, TRUE),
                    "VectorNoLog" = p4beta_(x, s1, s2, a, b),
                    "DoubleLog" = p4beta_(x[1], s1, s2, a, b, TRUE, TRUE),
                    "DoubleNoLog" = p4beta_(x[1], s1, s2, a, b),
                    "VectorLogNoLower" = p4beta_(x, s1, s2, a, b, FALSE, TRUE),
                    "VectorNoLogNoLower" = p4beta_(x, s1, s2, a, b, FALSE),
                    "DoubleLogNoLower" = p4beta_(x[1], s1, s2, a, b, FALSE, TRUE),
                    "DoubleNoLogNoLower" = p4beta_(x[1], s1, s2, a, b, FALSE)
                    )
                )
})


test_that('The quantile functions provide correct answers', {
    q4beta_ <- function(p, shape1, shape2, a, b,
                        lower_tail = TRUE, log_p = FALSE) {
        return((b - a) * qbeta(p, shape1, shape2, lower.tail = lower_tail,
               log.p = log_p) + a)
    }
    x <- c(0, 0.5, 1)
    s1 <- 2.0; s2 <- 2.0; a <- -2.5; b <- 2.5;
    expect_equal(test_q4beta_nolog(x, s1, s2, a, b),
                 list(
                    "VectorNoLog" = q4beta_(x, s1, s2, a, b),
                    "DoubleNoLog" = q4beta_(x[1], s1, s2, a, b),
                    "VectorNoLogNoLower" = q4beta_(x, s1, s2, a, b, FALSE),
                    "DoubleNoLogNoLower" = q4beta_(x[1], s1, s2, a, b, FALSE)
                    )
                )
    x <- c(-1, -2, -10)
    expect_equal(test_q4beta_log(x, s1, s2, a, b),
                 list(
                    "VectorLog" = q4beta_(x, s1, s2, a, b, TRUE, TRUE),
                    "DoubleLog" = q4beta_(x[1], s1, s2, a, b, TRUE, TRUE),
                    "VectorLogNoLower" = q4beta_(x, s1, s2, a, b, FALSE, TRUE),
                    "DoubleLogNoLower" = q4beta_(x[1], s1, s2, a, b, FALSE, TRUE)
                    )
                )
})

