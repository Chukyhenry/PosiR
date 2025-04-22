#' Fit OLS Model for a Given Subset of Predictors
#'
#' Fits a linear model for a specified subset of predictors using `lm.fit`.
#'
#' @param X_full Full design matrix (n x p_design), including intercept if added.
#' @param y Response vector (n).
#' @param q_indices Vector of column indices for the model.
#' @return A named vector of coefficients.
#' @keywords internal
fit_model_q <- function(X_full, y, q_indices) {
  if (length(q_indices) == 0) {
    stop("Model must contain at least one predictor/intercept.")
  }
  X_q <- X_full[, q_indices, drop = FALSE]

  # Use lm.fit for efficiency
  fit <- stats::lm.fit(X_q, y)

  # Handle rank deficiency
  if (fit$rank < length(q_indices)) {
    warning(paste("Model with columns", paste(colnames(X_q), collapse=", "), "is rank deficient."),
            call. = FALSE)
  }

  coeffs <- stats::setNames(fit$coefficients, colnames(X_q))
  return(coeffs)
}
