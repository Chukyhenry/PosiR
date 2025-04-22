#' Calculate Simultaneous Confidence Intervals using Bootstrap (Algorithm 1)
#'
#' Implements Algorithm 1 from Kuchibhotla et al. (2022) via bootstrap
#' for simultaneous confidence intervals over coefficients within a universe
#' of linear models.
#'
#' @param X Design matrix (numeric matrix, n x p). Should not include an intercept column
#'          if 'add_intercept = TRUE'. Column names are recommended.
#' @param y Response vector (numeric vector, length n).
#' @param Q_universe A list specifying the universe of models. Each element of the list
#'                   should be a vector of column indices (from the design matrix
#'                   including the intercept if `add_intercept=TRUE`) defining a model.
#'                   The list should be named, as names are used for model IDs.
#'                   Example: `list(model1 = c(1, 2), model2 = c(1, 3, 4))` assuming intercept is col 1.
#' @param alpha Significance level (e.g., 0.05 for 95% confidence).
#' @param B Number of bootstrap samples. Should be reasonably large (e.g., >= 1000).
#' @param add_intercept Logical. If TRUE, an intercept term is added as the first column
#'                      of the design matrix. Defaults to TRUE.
#' @param bootstrap_method Type of bootstrap. Currently only "pairs" (resampling rows of (X, y))
#'                         is implemented.
#' @param use_pbapply Logical. If TRUE and the `pbapply` package is installed, uses progress
#'                    bars from `pbapply`. Otherwise, uses base R `txtProgressBar`. Defaults to TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing:
#'         - `intervals`: A data frame with columns `model_id`, `coefficient_name`,
#'                      `estimate` (\( \widehat{\theta}_{q,j} \)), `lower`, `upper`,
#'                      `psi_hat_nqj` (\( \widehat{\Psi}_{n,q,j} \)), `se_nqj` (\( \sqrt{\widehat{\Psi}_{n,q,j}/n} \)).
#'         - `K_alpha`: The computed \( (1-\alpha) \) quantile (\( \widehat{K}_{\alpha} \)) of the max-t statistics.
#'         - `alpha`: Significance level used.
#'         - `B`: Number of bootstrap samples used.
#'         - `n_valid_T_star_b`: Number of bootstrap samples yielding a finite \( T^{*b} \) value.
#'         - `T_star_b`: Vector of the B bootstrap max-t statistics (may contain NAs).
#'         - `bootstrap_method`: Bootstrap method used.
#'         - `warnings`: A list of warnings generated during model fitting or calculations.
#'         - `valid_bootstrap_counts`: A nested list showing the number of valid bootstrap
#'                                     estimates used for each \( \widehat{\Psi}_{n,q,j} \).
#'
#' @export
#' @importFrom stats lm.fit quantile setNames sd var qnorm # Added var, qnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @suggests pbapply # Added suggests tag for optional dependency
#'
#' @examples
#' \dontrun{
#' # Install pbapply for better progress bars: install.packages("pbapply")
#'
#' set.seed(123)
#' n <- 100
#' p <- 4 # Number of non-intercept predictors
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("X", 1:p)
#' beta_true <- c(1, 0.5, 0, -0.5, 0) # Intercept, X1, X2, X3, X4
#' X_design_true <- cbind(1, X)
#' colnames(X_design_true) <- c("(Intercept)", colnames(X))
#' y <- X_design_true %*% beta_true + rnorm(n, sd = 1.5)
#'
#' # Define universe of models (using indices including intercept)
#' Q_models <- list(
#'   m_int_x1 = 1:2,          # Intercept + X1
#'   m_int_x1_x3 = c(1, 2, 4), # Intercept + X1 + X3
#'   m_all = 1:(p+1)          # Intercept + X1:X4
#' )
#'
#' # Run with B=100 for speed; use B=1000+ for actual analysis
#' # Use pbapply if installed (default), otherwise falls back to txtProgressBar
#' results <- simultaneous_ci(X, y, Q_universe = Q_models, alpha = 0.05, B = 100)
#' print(paste("K_alpha:", round(results$K_alpha, 3)))
#' print(results$intervals)
#'
#' # Compare with unadjusted intervals (using same Psi_hat)
#' z_alpha_2 <- stats::qnorm(1 - 0.05 / 2) # Use stats::qnorm explicitly if needed
#' unadj_intervals <- results$intervals
#' unadj_intervals$lower <- results$intervals$estimate - z_alpha_2 * results$intervals$se_nqj
#' unadj_intervals$upper <- results$intervals$estimate + z_alpha_2 * results$intervals$se_nqj
#' print("Unadjusted Intervals (for comparison):")
#' print(unadj_intervals[, c("model_id", "coefficient_name", "lower", "upper")])
#' }
simultaneous_ci <- function(X, y, Q_universe, alpha = 0.05, B = 1000,
                            add_intercept = TRUE, bootstrap_method = "pairs",
                            use_pbapply = TRUE, ...) {

  # --- Input Validation ---
  # (Keep existing validation checks for X, y, Q_universe, alpha, B, bootstrap_method)
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix.")
  if (!is.vector(y) || !is.numeric(y)) stop("y must be a numeric vector.")
  n <- nrow(X)
  p_orig <- ncol(X)
  if (n != length(y)) stop("Number of rows in X must match length of y.")
  if (!is.list(Q_universe) || length(Q_universe) == 0) stop("Q_universe must be a non-empty list.")
  if (is.null(names(Q_universe)) || any(names(Q_universe) == "")) stop("Q_universe must be a named list with non-empty names.")
  if (!all(sapply(Q_universe, is.numeric))) stop("Elements of Q_universe must be numeric vectors (indices).")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1.")
  if (!is.numeric(B) || B < 1) stop("B must be a positive integer.")
  if (B < 50) warning("Number of bootstrap samples B is small (< 50), results may be unstable.")
  if (bootstrap_method != "pairs") stop("Currently only 'pairs' bootstrap is supported.")

  # Check for pbapply if requested
  pbapply_available <- FALSE
  if (use_pbapply) {
    if (requireNamespace("pbapply", quietly = TRUE)) {
      pbapply_available <- TRUE
    } else {
      warning("`pbapply` package not found and `use_pbapply=TRUE`. Using base R `txtProgressBar` instead.", call. = FALSE)
    }
  }

  # --- Preprocessing ---
  # (Keep existing preprocessing for X_design, colnames, Q_universe_named)
  if (add_intercept) {
    X_design <- cbind(1, X)
    current_colnames <- if (is.null(colnames(X))) paste0("V", 1:p_orig) else colnames(X)
    colnames(X_design) <- c("(Intercept)", current_colnames)
  } else {
    X_design <- X
    if (is.null(colnames(X))) {
      colnames(X_design) <- paste0("V", 1:p_orig)
    }
  }
  p_design <- ncol(X_design)
  design_colnames <- colnames(X_design)

  Q_universe_named <- list()
  for (q_name in names(Q_universe)) {
    q_indices <- Q_universe[[q_name]]
    if (any(q_indices < 1) || any(q_indices > p_design)) {
      stop(paste("Indices in Q_universe element '", q_name, "' must be between 1 and ", p_design, ".", sep = ""))
    }
    if (any(duplicated(q_indices))) {
      warning(paste("Duplicate indices found in Q_universe element '", q_name, "'. Using unique indices.", sep = ""), call. = FALSE)
      q_indices <- unique(q_indices)
    }
    Q_universe_named[[q_name]] <- design_colnames[q_indices]
  }

  # --- Store Warnings ---
  fit_warnings <- list()
  record_warning <- function(w, context) {
    # Avoid recording the same warning repeatedly in loops if identical
    msg <- paste(context, w$message)
    if (!any(sapply(fit_warnings, function(x) x$full_message == msg))) {
      key <- paste(context, length(fit_warnings) + 1)
      fit_warnings[[key]] <<- list(message = w$message, context = context, full_message = msg)
    }
  }

  # --- Original Estimates ---
  original_estimates <- list()
  message("Fitting models on original data...")
  # Progress bar setup
  if (pbapply_available) {
    pb_orig <- pbapply::pboptions(type = "timer") # Use timer style
    on.exit(pbapply::pboptions(pb_orig), add = TRUE) # Restore options on exit
    iter_orig <- pbapply::pblapply(names(Q_universe_named), function(q_name) q_name) # Use pblapply for iteration tracking
  } else {
    pb_orig <- utils::txtProgressBar(min = 0, max = length(Q_universe_named), style = 3)
    iter_orig <- names(Q_universe_named)
  }

  for (i in seq_along(iter_orig)) {
    q_name <- iter_orig[[i]] # Get name regardless of iterator type
    q_col_names <- Q_universe_named[[q_name]]
    q_indices <- match(q_col_names, design_colnames)

    tryCatch({
      withCallingHandlers({
        # Assumes fit_model_q is available (moved to utils.R)
        coeffs <- fit_model_q(X_design, y, q_indices)
        original_estimates[[q_name]] <- stats::setNames(coeffs[q_col_names], q_col_names)
      }, warning = function(w) {
        record_warning(w, paste("Original fit:", q_name))
        invokeRestart("muffleWarning")
      })
    }, error = function(e) {
      stop(paste("Error fitting model '", q_name, "' on original data: ", e$message, sep = ""))
    })
    if (!pbapply_available) utils::setTxtProgressBar(pb_orig, i)
  }
  if (!pbapply_available) close(pb_orig)
  if (pbapply_available) pbapply::pblapply(1, function(x) NULL) # Finalize pbapply bar if needed

  # --- Bootstrap Loop ---
  # (Keep est_storage initialization)
  est_storage <- list()
  for(q_name in names(Q_universe_named)){
    est_storage[[q_name]] <- list()
    if(!is.null(original_estimates[[q_name]])){
      for(coeff_name in names(original_estimates[[q_name]])){
        est_storage[[q_name]][[coeff_name]] <- numeric(B)
      }
    }
  }

  sqrt_n <- sqrt(n)
  message(paste("Running", B, "bootstrap samples..."))
  # Progress bar setup for bootstrap loop (manual update)
  if (pbapply_available) {
    pb_boot <- pbapply::startpb(0, B)
    on.exit(pbapply::closepb(pb_boot), add = TRUE)
  } else {
    pb_boot <- utils::txtProgressBar(min = 0, max = B, style = 3)
  }

  for (b in 1:B) {
    # 1. Generate bootstrap sample
    indices_b <- sample(1:n, size = n, replace = TRUE)
    X_b <- X_design[indices_b, , drop = FALSE]
    y_b <- y[indices_b]

    # 2. Compute bootstrap estimators
    for (q_name in names(Q_universe_named)) {
      q_col_names <- Q_universe_named[[q_name]]
      q_indices <- match(q_col_names, design_colnames)

      tryCatch({
        withCallingHandlers({
          # Assumes fit_model_q is available
          coeffs_b <- fit_model_q(X_b, y_b, q_indices)
          orig_coeff_names <- names(original_estimates[[q_name]])
          if(!is.null(orig_coeff_names)){
            for(coeff_name in orig_coeff_names) {
              est_b <- coeffs_b[coeff_name]
              est_storage[[q_name]][[coeff_name]][b] <- ifelse(is.finite(est_b), est_b, NA_real_)
            }
          }
        }, warning = function(w) {
          record_warning(w, paste("Bootstrap fit:", q_name, "sample", b))
          invokeRestart("muffleWarning")
        })
      }, error = function(e) {
        warning(paste("Error fitting model '", q_name, "' in bootstrap sample ", b, ": ", e$message, ". Storing NAs.", sep = ""), call. = FALSE)
        orig_coeff_names <- names(original_estimates[[q_name]])
        if(!is.null(orig_coeff_names)){
          for(coeff_name in orig_coeff_names) {
            est_storage[[q_name]][[coeff_name]][b] <- NA_real_
          }
        }
      })
    } # End loop over q
    # Update progress bar
    if (pbapply_available) pbapply::setpb(pb_boot, b) else utils::setTxtProgressBar(pb_boot, b)
  } # End bootstrap loop (b)
  if (!pbapply_available) close(pb_boot)

  # --- Calculate Psi_hat_nqj ---
  message("Calculating bootstrap variance estimates (Psi_hat_nqj)...")
  psi_hat_nqj <- list()
  valid_bs_counts <- list()

  for (q_name in names(Q_universe_named)) {
    psi_hat_nqj[[q_name]] <- list()
    valid_bs_counts[[q_name]] <- list()
    q_coeffs_orig <- original_estimates[[q_name]]
    if(is.null(q_coeffs_orig)) next

    for(coeff_name in names(q_coeffs_orig)) {
      theta_qj <- q_coeffs_orig[[coeff_name]]
      theta_qj_star_b <- est_storage[[q_name]][[coeff_name]]

      valid_indices <- is.finite(theta_qj_star_b) & is.finite(theta_qj)
      n_valid_bs <- sum(valid_indices)
      valid_bs_counts[[q_name]][[coeff_name]] <- n_valid_bs

      if (n_valid_bs >= 2) {
        # *** FIX 1: Correct variance calculation ***
        # Calculate sum of squares of diffs, divide by (n_valid - 1)
        # No centering by mean(diffs) needed, as per Algorithm 1 structure.
        diffs <- sqrt_n * (theta_qj_star_b[valid_indices] - theta_qj)
        psi_hat_nqj[[q_name]][[coeff_name]] <- sum(diffs^2) / (n_valid_bs - 1)
      } else {
        psi_hat_nqj[[q_name]][[coeff_name]] <- NA_real_
        warning(paste("Could not estimate Psi_hat for '", coeff_name, "' in model '", q_name,
                      "' due to insufficient valid bootstrap samples (", n_valid_bs, ").", sep=""), call. = FALSE)
      }
    }
  }

  # --- Compute T*b and K_alpha ---
  message("Calculating bootstrap max-t statistics (T_star_b)...")
  T_star_b <- numeric(B)

  # Progress bar setup
  if (pbapply_available) {
    pb_T <- pbapply::startpb(0, B)
    on.exit(pbapply::closepb(pb_T), add = TRUE)
  } else {
    pb_T <- utils::txtProgressBar(min = 0, max = B, style = 3)
  }

  for(b in 1:B) {
    max_t_stat_b <- -Inf

    for(q_name in names(Q_universe_named)) {
      q_coeffs_orig <- original_estimates[[q_name]]
      if(is.null(q_coeffs_orig)) next

      for(coeff_name in names(q_coeffs_orig)) {
        theta_qj <- q_coeffs_orig[[coeff_name]]
        theta_qj_star_b <- est_storage[[q_name]][[coeff_name]][b]
        psi_hat <- psi_hat_nqj[[q_name]][[coeff_name]]

        # *** FIX 5: Stricter threshold for T*b stability ***
        # Check if all values needed are finite and psi_hat is sufficiently > 0
        if (is.finite(theta_qj) && is.finite(theta_qj_star_b) &&
            is.finite(psi_hat) && psi_hat > 1e-8) { # Stricter threshold
          t_stat_qjb <- abs(sqrt_n * (psi_hat^(-0.5)) * (theta_qj_star_b - theta_qj))
          if (is.finite(t_stat_qjb) && t_stat_qjb > max_t_stat_b) {
            max_t_stat_b <- t_stat_qjb
          }
        }
      }
    }
    T_star_b[b] <- ifelse(is.finite(max_t_stat_b), max_t_stat_b, NA_real_)
    # Update progress bar
    if (pbapply_available) pbapply::setpb(pb_T, b) else utils::setTxtProgressBar(pb_T, b)
  } # End loop for T*b
  if (!pbapply_available) close(pb_T)

  # (Keep validation checks for T_star_b and K_alpha calculation)
  valid_T_star_b <- T_star_b[is.finite(T_star_b)]
  n_valid_T_star_b <- length(valid_T_star_b)

  if (n_valid_T_star_b == 0) {
    stop("Could not compute any valid bootstrap T* statistics. Cannot determine K_alpha. Check model fits and variance estimates.")
  }
  if (n_valid_T_star_b < B * 0.5) {
    warning(paste("More than 50% of bootstrap T* statistics were non-finite (", B - n_valid_T_star_b, " / ", B,"). ",
                  "Results may be unreliable. Check model fits, variance estimates, and potential collinearity.", sep=""), call. = FALSE)
  }

  K_alpha <- stats::quantile(valid_T_star_b, probs = (1 - alpha), type = 8, na.rm = TRUE)


  # --- Construct Confidence Intervals ---
  message("Constructing confidence intervals...")
  results_list <- list()
  idx <- 1
  psi_threshold <- 1e-12 # Threshold for considering psi_hat as zero for CI construction

  for (q_name in names(Q_universe_named)) {
    q_coeffs_orig <- original_estimates[[q_name]]
    if(is.null(q_coeffs_orig)) next

    for(coeff_name in names(q_coeffs_orig)) {
      theta_qj <- q_coeffs_orig[[coeff_name]]
      psi_hat <- psi_hat_nqj[[q_name]][[coeff_name]]

      lower <- NA_real_
      upper <- NA_real_
      se_nqj <- NA_real_

      if (is.finite(theta_qj) && is.finite(psi_hat) && psi_hat >= 0) {
        # *** FIX 3: Add warning for near-zero psi_hat ***
        if (psi_hat > psi_threshold) {
          se_nqj <- sqrt(psi_hat / n)
          margin <- K_alpha * se_nqj
          lower <- theta_qj - margin
          upper <- theta_qj + margin
        } else {
          # Psi_hat is near zero, collapse interval to point estimate
          lower <- theta_qj
          upper <- theta_qj
          se_nqj <- 0
          warning(paste("Psi_hat for '", coeff_name, "' in model '", q_name,
                        "' is near zero (", format(psi_hat, digits=3, scientific=TRUE), # Format for clarity
                        "). Interval set to point estimate.", sep = ""),
                  call. = FALSE)
        }
      }

      results_list[[idx]] <- data.frame(
        model_id = q_name,
        coefficient_name = coeff_name,
        estimate = theta_qj,
        lower = lower,
        upper = upper,
        psi_hat_nqj = psi_hat,
        se_nqj = se_nqj,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  intervals_df <- do.call(rbind, results_list)
  # Order results for consistency
  if (nrow(intervals_df) > 0) {
    intervals_df <- intervals_df[order(intervals_df$model_id, intervals_df$coefficient_name), ]
  }
  rownames(intervals_df) <- NULL

  # --- Output ---
  # (Keep existing output structure)
  output <- list(
    intervals = intervals_df,
    K_alpha = as.numeric(K_alpha),
    alpha = alpha,
    B = B,
    n_valid_T_star_b = n_valid_T_star_b,
    T_star_b = T_star_b,
    bootstrap_method = bootstrap_method,
    warnings = fit_warnings, # Note: This now captures warnings from fits and near-zero psi
    valid_bootstrap_counts = valid_bs_counts
  )

  class(output) <- c("simultaneous_ci_result", "list")

  message("Done.")
  return(output)
}
