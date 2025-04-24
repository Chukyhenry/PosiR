#' Fit OLS model using lm.fit (Internal Helper)
#'
#' A simple wrapper around lm.fit for efficiency within the bootstrap loop.
#' Handles potential rank deficiency by returning NA for affected coefficients.
#' Intended for internal use within the PosiR package.
#'
#' @param X_full Full design matrix (n x p_design), including intercept if added.
#'               Assumed to have unique column names by the time it's used here.
#' @param y Response vector (n).
#' @param q_indices Vector of column indices (numeric) for the specific model q
#'                  to be fitted, relative to `X_full`.
#'
#' @return A named vector of coefficients corresponding to the columns specified
#'         by `q_indices`. Coefficients for linearly dependent columns identified
#'         by `lm.fit` will be `NA`. If `lm.fit` fails completely or loses coefficient
#'         names, returns a vector of NAs with the expected names.
#'
#' @keywords internal
#' @importFrom stats lm.fit setNames
#'
fit_model_q <- function(X_full, y, q_indices) {

  if (length(q_indices) == 0) {
    return(stats::setNames(numeric(0), character(0)))
  }
  if (any(q_indices < 1) || any(q_indices > ncol(X_full))) {
    stop("Invalid column indices passed to fit_model_q.", call. = FALSE)
  }

  q_colnames <- colnames(X_full)[q_indices]

  if (anyDuplicated(q_colnames)) {
    warning(paste("Duplicate column names detected within the selected model:",
                  paste(q_colnames[duplicated(q_colnames)], collapse=", "),
                  ". Results may be ambiguous."), call. = FALSE)
  }

  X_q <- X_full[, q_indices, drop = FALSE]

  fit <- tryCatch({
    stats::lm.fit(X_q, y)
  }, error = function(e) {
    warning(paste("lm.fit failed for model with columns:", paste(q_colnames, collapse=", "),
                  "; Error:", e$message, ". Returning NAs."), call. = FALSE)
    list(coefficients = stats::setNames(rep(NA_real_, length(q_colnames)), q_colnames))
  })

  coeffs <- fit$coefficients

  output_coeffs <- stats::setNames(rep(NA_real_, length(q_colnames)), q_colnames)

  if (!is.null(names(coeffs))) {
    valid_fit_coeffs <- coeffs[!is.na(coeffs)]
    common_names <- intersect(names(output_coeffs), names(valid_fit_coeffs))
    if (length(common_names) > 0) {
      output_coeffs[common_names] <- valid_fit_coeffs[common_names]
    }
    missing_names <- setdiff(names(output_coeffs), names(coeffs))
    if (length(missing_names) > 0) {
      warning(paste("Coefficients missing from lm.fit result:", paste(missing_names, collapse=", ")), call.=FALSE)
    }

  } else if (length(coeffs) > 0) {
    warning("Coefficient names were lost during fitting. Cannot reliably assign results. Returning NAs.", call.=FALSE)
  }

  return(output_coeffs)
}
