context("Location-Scale t Distribution")

test_that('The density functions provide correct answers', {
    dlst_ <- function(x, df, mu, sigma, log_p = FALSE) {
        p = dt((x - mu)/sigma, df, log = log_p)
        if ( log_p ) {
            return(p - log(sigma))
        }
        return(p / sigma)
    }
    x <- c(-3, 2, 0, 4, -1)
    df <- 1; mu <- 1; sigma <- 2;
    expect_equal(test_dlst(x, df, mu, sigma),
                 list(
                    "VectorLog" = dlst_(x, df, mu, sigma, TRUE),
                    "VectorNoLog" = dlst_(x, df, mu, sigma, FALSE),
                    "DoubleLog" = dlst_(x[1], df, mu, sigma, TRUE),
                    "DoubleNoLog" = dlst_(x[1], df, mu, sigma, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_dlst(x, df, mu, sigma),
                 list(
                    "VectorLog" = dlst_(x, df, mu, sigma, TRUE),
                    "VectorNoLog" = dlst_(x, df, mu, sigma, FALSE),
                    "DoubleLog" = dlst_(x[1], df, mu, sigma, TRUE),
                    "DoubleNoLog" = dlst_(x[1], df, mu, sigma, FALSE)
                    )
                )
})


test_that('The distribution functions provide correct answers', {
    plst_ <- function(q, df, mu, sigma, lower_tail = TRUE, log_p = FALSE) {
        return(pt((q - mu)/sigma, df, lower.tail = lower_tail, log.p = log_p))
    }
    x <- c(-3, 2, 0, 4, -1)
    df <- 1; mu <- 1; sigma <- 2;
    expect_equal(test_plst(x, df, mu, sigma),
                 list(
                    "VectorLog" = plst_(x, df, mu, sigma, TRUE, TRUE),
                    "VectorNoLog" = plst_(x, df, mu, sigma),
                    "DoubleLog" = plst_(x[1], df, mu, sigma, TRUE, TRUE),
                    "DoubleNoLog" = plst_(x[1], df, mu, sigma),
                    "VectorLogNoLower" = plst_(x, df, mu, sigma, FALSE, TRUE),
                    "VectorNoLogNoLower" = plst_(x, df, mu, sigma, FALSE),
                    "DoubleLogNoLower" = plst_(x[1], df, mu, sigma, FALSE, TRUE),
                    "DoubleNoLogNoLower" = plst_(x[1], df, mu, sigma, FALSE)
                    )
                )
    x <- x[-1]
    expect_equal(test_plst(x, df, mu, sigma),
                 list(
                    "VectorLog" = plst_(x, df, mu, sigma, TRUE, TRUE),
                    "VectorNoLog" = plst_(x, df, mu, sigma),
                    "DoubleLog" = plst_(x[1], df, mu, sigma, TRUE, TRUE),
                    "DoubleNoLog" = plst_(x[1], df, mu, sigma),
                    "VectorLogNoLower" = plst_(x, df, mu, sigma, FALSE, TRUE),
                    "VectorNoLogNoLower" = plst_(x, df, mu, sigma, FALSE),
                    "DoubleLogNoLower" = plst_(x[1], df, mu, sigma, FALSE, TRUE),
                    "DoubleNoLogNoLower" = plst_(x[1], df, mu, sigma, FALSE)
                    )
                )
})


test_that('The quantile functions provide correct answers', {
    qlst_ <- function(p, df, mu, sigma, lower_tail = TRUE, log_p = FALSE) {
        return(qt(p, df, lower.tail = lower_tail, log.p = log_p) * sigma + mu)
    }
    x <- c(0, 0.5, 1)
    df <- 1; mu <- 1; sigma <- 2;
    expect_equal(test_qlst_nolog(x, df, mu, sigma),
                 list(
                    "VectorNoLog" = qlst_(x, df, mu, sigma),
                    "DoubleNoLog" = qlst_(x[1], df, mu, sigma),
                    "VectorNoLogNoLower" = qlst_(x, df, mu, sigma, FALSE),
                    "DoubleNoLogNoLower" = qlst_(x[1], df, mu, sigma, FALSE)
                    )
                )
    x <- c(-1, -2, -10)
    expect_equal(test_qlst_log(x, df, mu, sigma),
                 list(
                    "VectorLog" = qlst_(x, df, mu, sigma, TRUE, TRUE),
                    "DoubleLog" = qlst_(x[1], df, mu, sigma, TRUE, TRUE),
                    "VectorLogNoLower" = qlst_(x, df, mu, sigma, FALSE, TRUE),
                    "DoubleLogNoLower" = qlst_(x[1], df, mu, sigma, FALSE, TRUE)
                    )
                )
})

