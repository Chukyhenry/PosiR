
test_that("Basic functionality works with intercept-only model", {
  X <- matrix(1, 10, 1, dimnames = list(NULL, "(Intercept)"))
  y <- rnorm(10)
  Q <- list(null_model = 1)
  result <- simultaneous_ci(X, y, Q, B = 50, add_intercept = FALSE, verbose = FALSE)
  expect_s3_class(result, "simultaneous_ci_result")
  expect_true(nrow(result$intervals) == 1)
})

test_that("Function adds intercept and handles single predictor", {
  X <- matrix(rnorm(10), 10, 1)
  colnames(X) <- "X1"
  y <- rnorm(10)
  Q <- list(model_with_intercept = 1:2)  # intercept will be col 1, X1 will be col 2
  result <- simultaneous_ci(X, y, Q, B = 50, add_intercept = TRUE, verbose = FALSE)
  expect_s3_class(result, "simultaneous_ci_result")
  expect_true("X1" %in% result$intervals$coefficient_name)
})

test_that("Handles perfect collinearity gracefully", {
  set.seed(123)
  x1 <- rnorm(20)
  x2 <- x1  # perfectly collinear
  X <- cbind(1, x1, x2)
  colnames(X) <- c("(Intercept)", "V1", "V2")

  y <- rnorm(20)
  Q <- list(collinear = 1:3)

  result <- simultaneous_ci(X, y, Q, B = 50, add_intercept = FALSE, verbose = FALSE)

  # Check: at least one variance estimate is NA or very small
  psi_values <- unlist(result$intervals$psi_hat_nqj)
  expect_true(any(is.na(psi_values) | psi_values < 1e-6))
})


test_that("Returns proper structure with valid input", {
  X <- matrix(rnorm(100*2), 100, 2)
  colnames(X) <- c("X1", "X2")  # <- add unique names
  y <- X[,1] + rnorm(100)
  Q <- list(m1 = 1:2)
  result <- simultaneous_ci(X, y, Q, B = 100, verbose = FALSE)
  expect_named(result,
               c("intervals", "K_alpha", "alpha", "B", "n_valid_T_star_b",
                 "T_star_b", "bootstrap_method", "warnings_list",
                 "valid_bootstrap_counts", "n_bootstrap_errors")
  )
})

test_that("Handles multiple models in Q_universe", {
  X <- matrix(rnorm(100*3), 100, 3, dimnames = list(NULL, c("X1", "X2", "X3")))
  y <- rnorm(100)
  Q <- list(
    model1 = 1:2,           # (Intercept) + X1
    model2 = 1:3,           # (Intercept) + X1 + X2
    model3 = c(1, 4)        # (Intercept) + X3
  )
  result <- simultaneous_ci(X, y, Q, B = 100, add_intercept = TRUE, verbose = FALSE)

  # Verify all models appear in results
  expect_equal(length(unique(result$intervals$model_id)), 3)

  # Verify coefficient counts
  expect_equal(nrow(result$intervals),
               2 + 3 + 2) # 2 coefficients in model1, 3 in model2, 2 in model3
})

test_that("Sequential mode works and produces output", {
  X <- matrix(rnorm(200*2), 200, 2, dimnames = list(NULL, c("V1", "V2")))
  y <- X[,1] + rnorm(200)
  Q <- list(full_model = 1:3) # Assumes add_intercept=TRUE

  set.seed(123)
  res_seq <- simultaneous_ci(X, y, Q, B = 100, cores = 1, seed = 123)

  expect_s3_class(res_seq, "simultaneous_ci_result")
})

test_that("Parallel mode works and produces consistent results", {
  X <- matrix(rnorm(200*2), 200, 2, dimnames = list(NULL, c("V1", "V2")))
  y <- X[,1] + rnorm(200)
  Q <- list(full_model = 1:3)

  set.seed(123)
  res_par <- simultaneous_ci(X, y, Q, B = 100, cores = 2, seed = 123)

  expect_s3_class(res_par, "simultaneous_ci_result")
})

test_that("Parallel and sequential modes are approximately consistent", {
  X <- matrix(rnorm(200*2), 200, 2, dimnames = list(NULL, c("V1", "V2")))
  y <- X[,1] + rnorm(200)
  Q <- list(full_model = 1:3)

  set.seed(123)
  res_seq <- simultaneous_ci(X, y, Q, B = 100, cores = 1, seed = 123)

  set.seed(123)
  res_par <- simultaneous_ci(X, y, Q, B = 100, cores = 2, seed = 123)

  # Only compare critical values (allow slight randomness tolerance)
  expect_equal(res_seq$K_alpha, res_par$K_alpha, tolerance = 0.25)
})

test_that("Handles small sample sizes with warnings", {
  X <- matrix(rnorm(5*2), 5, 2, dimnames = list(NULL, c("X1", "X2")))
  y <- rnorm(5)
  Q <- list(small_model = 1:2)

  expect_warning(
    simultaneous_ci(X, y, Q, B = 50, verbose = FALSE),
    "Sample size.*unstable"
  )
})

test_that("Rejects invalid bootstrap methods", {
  X <- matrix(rnorm(100*2), 100, 2, dimnames = list(NULL, c("X1", "X2")))
  y <- rnorm(100)
  Q <- list(m1 = 1:2)

  expect_error(
    simultaneous_ci(X, y, Q, bootstrap_method = "invalid"),
    "only.*pairs.*bootstrap", ignore.case = TRUE
  )
})


test_that("Detects duplicate model names", {
  X <- matrix(rnorm(100*2), 100, 2, dimnames = list(NULL, c("X1", "X2")))
  y <- rnorm(100)
  Q <- list(1:2, 1:3)
  names(Q) <- c("m1", "m1")  # Force duplicate names

  expect_error(
    simultaneous_ci(X, y, Q),
    "Duplicate model names detected in Q_universe"
  )
})

test_that("Bootstrap variances approximate theoretical values", {
  # Simple 1-predictor scenario with known variance
  X <- matrix(1:100, ncol = 1)
  colnames(X) <- "X1"
  y <- 2*X[,1] + rnorm(100, sd = 5)
  Q <- list(simple_model = 1:2) # Intercept + X1

  result <- simultaneous_ci(X, y, Q, B = 1000, verbose = FALSE)

  # Theoretical SE for slope = sigma/sqrt(Sxx)
  sigma_hat <- sd(y - 2*X[,1])
  Sxx <- sum((X[,1] - mean(X[,1]))^2)
  theoretical_se <- sigma_hat/sqrt(Sxx)

  # Compare with bootstrap SE (allow 20% tolerance)
  bootstrap_se <- result$intervals$se_nqj[2] # Slope term
  expect_equal(bootstrap_se, theoretical_se, tolerance = 0.2)
})

test_that("Simultaneous intervals achieve nominal coverage", {
  true_beta <- 0.5
  n_sim <- 50
  B <- 200
  alpha <- 0.05
  coverage_tolerance <- 0.05  # Increased tolerance for small n_sim

  coverage <- replicate(n_sim, {
    X <- matrix(rnorm(100*2), 100, 2, dimnames = list(NULL, c("V1", "V2")))
    y <- X[,1] * true_beta + rnorm(100)
    Q <- list(true_model = 1:2) # Intercept + V1 + V2

    result <- simultaneous_ci(X, y, Q, B = B, alpha = alpha, verbose = FALSE)

    # Check coverage for ALL parameters simultaneously
    all_covered <- TRUE
    for (coef_row in 1:nrow(result$intervals)) {
      coef_name <- result$intervals$coefficient_name[coef_row]
      true_value <- ifelse(coef_name == "V1", true_beta, 0) # Intercept truth = 0
      interval <- result$intervals[coef_row, c("lower", "upper")]
      all_covered <- all_covered && (interval$lower <= true_value && true_value <= interval$upper)
    }
    all_covered
  })

  coverage_rate <- mean(coverage)
  message("Simultaneous coverage rate: ", round(coverage_rate, 3))
  expect_true(coverage_rate >= (1 - alpha) - coverage_tolerance) # Allow >=90% coverage
}) |> suppressWarnings()

test_that("Marginal coverage matches expectations", {
  true_beta <- 0.5
  n_sim <- 100
  B <- 500
  alpha <- 0.05
  coverage_tolerance <- 0.05

  coverage <- replicate(n_sim, {
    X <- matrix(rnorm(100), 100, 1, dimnames = list(NULL, "V1"))
    y <- X[,1] * true_beta + rnorm(100)
    Q <- list(single_param = 1) # Only V1, no intercept

    result <- simultaneous_ci(X, y, Q, B = B, add_intercept = FALSE, verbose = FALSE)
    interval <- result$intervals[result$intervals$coefficient_name == "V1", c("lower", "upper")]
    interval$lower <= true_beta && true_beta <= interval$upper
  })

  coverage_rate <- mean(coverage)
  message("Marginal coverage rate: ", round(coverage_rate, 3))
  expect_true(abs(coverage_rate - (1 - alpha)) < coverage_tolerance)
})

test_that("Plot method returns invisibly", {
  X <- matrix(rnorm(100*2), 100, 2, dimnames = list(NULL, c("Age", "Income")))
  y <- rnorm(100)
  Q <- list(model = 1:2)
  result <- simultaneous_ci(X, y, Q, B = 50, verbose = FALSE)

  # Only test invisibility here, do not assign
  expect_invisible(plot(result, verbose = FALSE))
})
