# tests/testthat/test-simultaneous_ci.R
test_that("Basic functionality works with intercept-only model", {
  X <- matrix(1, 10, 1, dimnames = list(NULL, "(Intercept)"))
  y <- rnorm(10)
  Q <- list(null_model = 1)
  result <- simultaneous_ci(X, y, Q, B = 50, verbose = FALSE)
  expect_s3_class(result, "simultaneous_ci_result")
  expect_true(nrow(result$intervals) == 1)
})

test_that("Handles perfect collinearity gracefully", {
  X <- cbind(1, MASS::mvrnorm(20, mu=c(0,0), Sigma=matrix(c(1,1,1,1), ncol=2)),
             colnames(X) <- c("(Intercept)", "V1", "V2"),
             y <- rnorm(20),
             Q <- list(collinear = 1:3),
             expect_warning(
               simultaneous_ci(X, y, Q, B = 50, verbose = FALSE),
               regexp = "insufficient valid bootstrap samples|rank-deficient")
)
})
  test_that("Returns proper structure with valid input", {
    X <- matrix(rnorm(100*2), 100, 2)
    y <- X[,1] + rnorm(100)
    Q <- list(m1 = 1:2)
    result <- simultaneous_ci(X, y, Q, B = 100, verbose = FALSE)
    expect_named(result,
                 c("intervals", "K_alpha", "alpha", "B", "n_valid_T_star_b",
                   "T_star_b", "bootstrap_method", "warnings_list",
                   "valid_bootstrap_counts", "n_bootstrap_errors")
    )
  })
