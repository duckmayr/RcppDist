context("Inverse Wishart Distribution")

test_that('The density function provides correct answers', {
    S = matrix(c(1/2, 1/3, 1/3, 1/4), nrow = 2)
    expect_equal(test_diwish(S, 2, S),
                 list(
                    "Log" = 2.88397493155479,
                    "NoLog" = 17.885224616332
                    )
                )
})
