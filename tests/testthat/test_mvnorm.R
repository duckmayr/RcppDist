context("Multivariate Normal Distribution")

test_that('The density function provides correct answers', {
    dmvnorm_ <- function(x, mu, S, log = FALSE) {
        n = nrow(x)
        det_S = det(S)
        S = solve(S)
        result = numeric(n)
        if ( log ) {
            P = -1.0 * (ncol(x) / 2.0) * log(2 * pi) - 0.5 * log(det_S)
            for ( i in 1:n ) {
                X = x[i, ] - mu
                result[i] = P - 0.5 * t(X) %*% S %*% X
            }
        }
        else {
            P = 1.0 / sqrt(((2 * pi)^nrow(S)) * det_S)
            for ( i in 1:n ) {
                X = x[i, ] - mu
                result[i] = P * exp(-0.5 * t(X) %*% S %*% X)
            }
        }
        return(result)
    }
    mu = c(1, -1); S = matrix(c(1/2, 1/3, 1/3, 1/4), nrow = 2)
    x = matrix(runif(20, -5, 5), nrow = 10, ncol = 2)
    expect_equal(test_dmvnorm(x, mu, S),
                 list(
                    "Log" = matrix(dmvnorm_(x, mu, S, TRUE), ncol = 1),
                    "NoLog" = matrix(dmvnorm_(x, mu, S, FALSE), ncol = 1)
                    )
                )
})
