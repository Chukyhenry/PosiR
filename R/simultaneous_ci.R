#' Compute Simultaneous Confidence Intervals via Bootstrap (Post-Selection Inference)
#'
#' Implements Algorithm 1 from Kuchibhotla et al. (2022) using bootstrap-based max-t statistics
#' to construct valid simultaneous confidence intervals for selected regression coefficients
#' across a user-specified universe of linear models.
#'
#' Supports parallel execution, internal warnings capture, and returns structured results
#' with estimates, intervals, bootstrap diagnostics, and inference statistics.
#'
#' @param X Numeric matrix (n x p): Design matrix. Must have unique column names.
#'          Do not include an intercept if `add_intercept = TRUE`.
#' @param y Numeric vector (length n): Response vector.
#' @param Q_universe Named list of numeric vectors. Each element specifies a model as a
#'        vector of column indices (accounting for intercept if `add_intercept = TRUE`).
#'        Names are used to identify each model in results.
#' @param alpha Significance level for the confidence intervals. Default is 0.05.
#' @param B Integer. Number of bootstrap samples. Default is 1000.
#' @param add_intercept Logical. If TRUE, adds an intercept as the first column of the design matrix. Default is TRUE.
#' @param bootstrap_method Character. Bootstrap type. Only "pairs" is currently supported.
#' @param cores Integer. Number of CPU cores to use for bootstrap parallelization. Default is 1.
#' @param use_pbapply Logical. Use `pbapply` for progress bars if available. Default is TRUE.
#' @param seed Optional numeric. Random seed for reproducibility. Used for parallel-safe RNG.
#' @param verbose Logical. Whether to display status messages. Default is TRUE.
#' @param ... Reserved for future use.
#'
#' @return A list of class `simultaneous_ci_result` with elements:
#'
#' - `intervals`: Data frame with estimates, confidence intervals, variances, and SEs
#' - `K_alpha`: Bootstrap (1 - alpha) quantile of max-t statistics
#' - `T_star_b`: Vector of bootstrap max-t statistics
#' - `n_valid_T_star_b`: Number of finite bootstrap max-t statistics
#' - `alpha`, `B`, `bootstrap_method`: Metadata
#' - `warnings_list`: Internal warnings collected during bootstrap/model fitting
#' - `valid_bootstrap_counts`: Valid bootstrap replicates per parameter
#' - `n_bootstrap_errors`: Total bootstrap fitting errors
#'
#' @references
#' Kuchibhotla, A., Kolassa, J., & Kuffner, T. (2022). Post-selection inference.
#' *Annual Review of Statistics and Its Application*, 9(1), 505–527.
#'
#' @export
#' @importFrom stats lm.fit quantile setNames sd var qnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @import parallel
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(100 * 2), 100, 2, dimnames = list(NULL, c("X1", "X2")))
#' y <- X[,1] * 0.5 + rnorm(100)
#' Q <- list(model = 1:2)
#' res <- simultaneous_ci(X, y, Q, B = 100, cores = 1)
#' print(res$intervals)
#' plot(res)
#' }
simultaneous_ci <- function(X, y, Q_universe, alpha = 0.05, B = 1000,
                            add_intercept = TRUE, bootstrap_method = "pairs",
                            cores = 1, use_pbapply = TRUE, seed = NULL, verbose = TRUE, ...) {

  # --- Input Validation ---
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix.")
  if (is.null(colnames(X)) || anyDuplicated(colnames(X))) {
    stop("Input matrix 'X' must have unique column names.")
  }
  if (!is.vector(y) || !is.numeric(y)) stop("y must be a numeric vector.")
  n <- nrow(X)
  p_orig <- ncol(X)
  if (n != length(y)) stop("Number of rows in X must match length of y.")
  if (!is.list(Q_universe) || length(Q_universe) == 0) stop("Q_universe must be a non-empty list.")
  if (is.null(names(Q_universe)) || any(names(Q_universe) == "")) {
    stop("Q_universe must be a named list with non-empty names.")
  }
  if (anyDuplicated(names(Q_universe))) {
    stop("Duplicate model names detected in Q_universe.")
  }
  if (!all(sapply(Q_universe, is.numeric))) stop("Elements of Q_universe must be numeric vectors (indices).")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1.")
  if (!is.numeric(B) || B < 1) stop("B must be a positive integer.")
  if (B < 50) warning("Number of bootstrap samples B is small (< 50), results may be unstable.")
  if (bootstrap_method != "pairs") stop("Currently only 'pairs' bootstrap is supported.")
  if (!is.logical(verbose) || length(verbose) != 1) stop("verbose must be TRUE or FALSE.")

  if (!is.numeric(cores) || cores < 1) stop("cores must be a positive integer.")
  cores <- as.integer(cores)
  available_cores <- parallel::detectCores()
  if (cores > available_cores) {
    warning(paste("cores requested (", cores, ") exceeds available cores (", available_cores, "). Using ", available_cores, " cores.", sep=""), call. = FALSE)
    cores <- available_cores
  }

  pbapply_available <- FALSE
  if (use_pbapply) {
    if (requireNamespace("pbapply", quietly = TRUE)) {
      pbapply_available <- TRUE
    } else {
      if (verbose) warning("`pbapply` package not found and `use_pbapply=TRUE`. Progress bars disabled for parallel or using base R for sequential.", call. = FALSE)
    }
  }

  # --- Preprocessing ---
  if (add_intercept) {
    X_design <- cbind(1, X)
    colnames(X_design) <- c("(Intercept)", colnames(X))
  } else {
    X_design <- X
  }
  p_design <- ncol(X_design)
  design_colnames <- colnames(X_design)
  n <- nrow(X_design)
  if (n < 10) {
    warning("Sample size is very small (n < 10); results may be unstable.", call. = FALSE)
  }

  if (add_intercept && "(Intercept)" %in% colnames(X)) {
    warning("Input matrix 'X' already contains a column named '(Intercept)' and 'add_intercept=TRUE'. Ensure this is intended.", call. = FALSE)
  }
  if (anyDuplicated(design_colnames)) {
    stop("Internal Error: Column names in design matrix (X_design) are not unique after adding intercept. Check input 'X' names.")
  }

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
  warnings_list <- list()
  add_warning <- function(msg, context = "", call = FALSE) {
    warnings_list[[length(warnings_list) + 1]] <<- list(
      context = context,
      message = msg
    )
    warning(paste0("[", context, "] ", msg), call. = call)
  }

  # --- Original Estimates ---
  original_estimates <- list()
  if (verbose) message("Fitting models on original data...")
  for (q_name in names(Q_universe_named)) {
    q_col_names <- Q_universe_named[[q_name]]
    q_indices <- match(q_col_names, design_colnames)
    tryCatch({
      withCallingHandlers({
        coeffs <- fit_model_q(X_design, y, q_indices)
        original_estimates[[q_name]] <- stats::setNames(coeffs[q_col_names], q_col_names)
      }, warning = function(w) {
        add_warning(w$message, paste("Original fit:", q_name))
        invokeRestart("muffleWarning")
      })
    }, error = function(e) {
      stop(paste("Error fitting model '", q_name, "' on original data: ", e$message, sep = ""))
    })
  }

  # --- Bootstrap Loop ---
  sqrt_n <- sqrt(n)
  if (verbose) message(paste("Running", B, "bootstrap samples", ifelse(cores > 1, paste("on", cores, "cores"), ""), "..."))

  bootstrap_iteration <- function(b_index, .X_design, .y, .n, .Q_universe_named, .design_colnames) {
    iter_results <- list(estimates = list(), error_occurred = FALSE, warnings = list())
    local_warnings <- list()
    local_add_warning <- function(msg, context){
      local_warnings[[length(local_warnings) + 1]] <<- list(message=msg, context=context)
    }
    tryCatch({
      indices_b <- sample(1:.n, size = .n, replace = TRUE)
      X_b <- .X_design[indices_b, , drop = FALSE]
      y_b <- .y[indices_b]
      iter_estimates <- list()
      for (q_name in names(.Q_universe_named)) {
        q_col_names <- .Q_universe_named[[q_name]]
        q_indices <- match(q_col_names, .design_colnames)
        withCallingHandlers({
          coeffs_b <- fit_model_q(X_b, y_b, q_indices)
          iter_estimates[[q_name]] <- stats::setNames(coeffs_b[q_col_names], q_col_names)
        }, warning = function(w) {
          local_add_warning(w$message, paste("Bootstrap fit:", q_name, "sample", b_index))
          invokeRestart("muffleWarning")
        })
      }
      iter_results$estimates <- iter_estimates
    }, error = function(e) {
      iter_results$error_occurred <- TRUE
      local_add_warning(e$message, paste("Bootstrap fit error sample", b_index))
    })
    iter_results$warnings <- local_warnings
    return(iter_results)
  }

  bootstrap_results <- list()
  cl <- NULL

  if (cores > 1) {
    if (verbose) message("Setting up parallel cluster...")
    cl <- parallel::makeCluster(cores)
    on.exit({
      if (!is.null(cl)) parallel::stopCluster(cl)
    }, add = TRUE)  # << safer stopCluster only if cl exists

    if (!is.null(seed) && is.numeric(seed)) {
      parallel::clusterSetRNGStream(cl, seed)
    } else {
      warning("No seed provided for parallel execution. Results may not be exactly reproducible.", call. = FALSE)
    }

    parallel::clusterExport(cl, varlist = c("fit_model_q"), envir = environment(fit_model_q))

    if (pbapply_available && verbose) {
      bootstrap_results <- pbapply::pblapply(1:B, bootstrap_iteration,
                                             .X_design=X_design, .y=y, .n=n,
                                             .Q_universe_named=Q_universe_named,
                                             .design_colnames=design_colnames, cl=cl)
    } else {
      bootstrap_results <- parallel::parLapply(cl, 1:B, bootstrap_iteration,
                                               .X_design=X_design, .y=y, .n=n,
                                               .Q_universe_named=Q_universe_named,
                                               .design_colnames=design_colnames)
    }
    parallel::stopCluster(cl); cl <- NULL

  } else {
    if (!is.null(seed) && is.numeric(seed)) { set.seed(seed) }
    if (pbapply_available && verbose) {
      bootstrap_results <- pbapply::pblapply(1:B, bootstrap_iteration,
                                             .X_design=X_design, .y=y, .n=n,
                                             .Q_universe_named=Q_universe_named,
                                             .design_colnames=design_colnames)
    } else {
      pb_boot <- NULL
      if (verbose) pb_boot <- utils::txtProgressBar(min = 0, max = B, style = 3)
      results_temp <- vector("list", B)
      for(b in 1:B){
        results_temp[[b]] <- bootstrap_iteration(b, .X_design=X_design, .y=y, .n=n,
                                                 .Q_universe_named=Q_universe_named,
                                                 .design_colnames=design_colnames)
        if (verbose && !is.null(pb_boot)) utils::setTxtProgressBar(pb_boot, b)
      }
      if (verbose && !is.null(pb_boot)) close(pb_boot)
      bootstrap_results <- results_temp
    }
  }

  # --- Process Bootstrap Results ---
  if (verbose) message("Processing bootstrap results...")
  total_bootstrap_errors <- 0
  for (res in bootstrap_results) {
    if (res$error_occurred) total_bootstrap_errors <- total_bootstrap_errors + 1
    for(w in res$warnings){
      add_warning(w$message, w$context)
    }
  }
  if (total_bootstrap_errors > 0) {
    add_warning(paste(total_bootstrap_errors, "out of", B, "bootstrap iterations encountered an error during model fitting."), "Bootstrap Summary")
  }

  est_storage <- list()
  for(q_name in names(Q_universe_named)){
    est_storage[[q_name]] <- list()
    if(!is.null(original_estimates[[q_name]])){
      for(coeff_name in names(original_estimates[[q_name]])){
        est_storage[[q_name]][[coeff_name]] <- rep(NA_real_, B)
      }
    }
  }
  for (b in 1:B) {
    iter_estimates <- bootstrap_results[[b]]$estimates
    if (bootstrap_results[[b]]$error_occurred || is.null(iter_estimates)) next
    for (q_name in names(iter_estimates)) {
      if (!is.null(est_storage[[q_name]])){
        q_iter_coeffs <- iter_estimates[[q_name]]
        for (coeff_name in names(q_iter_coeffs)) {
          if (!is.null(est_storage[[q_name]][[coeff_name]])){
            est_b <- q_iter_coeffs[[coeff_name]]
            est_storage[[q_name]][[coeff_name]][b] <- ifelse(is.finite(est_b), est_b, NA_real_)
          }
        }
      }
    }
  }

  # --- Calculate Psi_hat_nqj ---
  if (verbose) message("Calculating bootstrap variance estimates (Psi_hat_nqj)...")
  psi_hat_nqj <- list()
  valid_bs_counts <- list()
  for (q_name in names(Q_universe_named)) {
    psi_hat_nqj[[q_name]] <- list()
    valid_bs_counts[[q_name]] <- list()
    q_coeffs_orig <- original_estimates[[q_name]]
    if(is.null(q_coeffs_orig)) next
    for(coeff_name in names(q_coeffs_orig)) {
      theta_qj <- q_coeffs_orig[[coeff_name]]
      if(is.null(est_storage[[q_name]]) || is.null(est_storage[[q_name]][[coeff_name]])){
        psi_hat_nqj[[q_name]][[coeff_name]] <- NA_real_
        valid_bs_counts[[q_name]][[coeff_name]] <- 0
        next
      }
      theta_qj_star_b <- est_storage[[q_name]][[coeff_name]]
      valid_indices <- is.finite(theta_qj_star_b) & is.finite(theta_qj)
      n_valid_bs <- sum(valid_indices)
      valid_bs_counts[[q_name]][[coeff_name]] <- n_valid_bs
      if (n_valid_bs >= 2) {
        diffs <- sqrt_n * (theta_qj_star_b[valid_indices] - theta_qj)
        psi_hat_nqj[[q_name]][[coeff_name]] <- sum(diffs^2) / (n_valid_bs - 1)
      } else {
        psi_hat_nqj[[q_name]][[coeff_name]] <- NA_real_
        if(is.finite(theta_qj)){
          add_warning(paste("Could not estimate Psi_hat for '", coeff_name, "' in model '", q_name,
                            "' due to insufficient valid bootstrap samples (", n_valid_bs, ").", sep=""), "Psi Calculation")
        }
      }
    }
  }

  # --- Compute T*b and K_alpha ---
  if (verbose) message("Calculating bootstrap max-t statistics (T_star_b)...")
  T_star_b <- numeric(B)
  psi_t_stat_threshold <- 1e-8
  for(b in 1:B) {
    max_t_stat_b <- -Inf
    if (bootstrap_results[[b]]$error_occurred) {
      T_star_b[b] <- NA_real_
      next
    }
    iter_estimates <- bootstrap_results[[b]]$estimates
    for(q_name in names(iter_estimates)) {
      q_coeffs_orig <- original_estimates[[q_name]]
      if(is.null(q_coeffs_orig)) next
      q_iter_coeffs <- iter_estimates[[q_name]]
      for(coeff_name in names(q_iter_coeffs)) {
        if(is.null(q_coeffs_orig[[coeff_name]]) || is.null(psi_hat_nqj[[q_name]][[coeff_name]])) next
        theta_qj <- q_coeffs_orig[[coeff_name]]
        theta_qj_star_b <- q_iter_coeffs[[coeff_name]]
        psi_hat <- psi_hat_nqj[[q_name]][[coeff_name]]
        if (is.finite(theta_qj) && is.finite(theta_qj_star_b) &&
            is.finite(psi_hat) && psi_hat > psi_t_stat_threshold) {
          t_stat_qjb <- abs(sqrt_n * (psi_hat^(-0.5)) * (theta_qj_star_b - theta_qj))
          if (is.finite(t_stat_qjb) && t_stat_qjb > max_t_stat_b) {
            max_t_stat_b <- t_stat_qjb
          }
        }
      }
    }
    T_star_b[b] <- ifelse(is.finite(max_t_stat_b), max_t_stat_b, NA_real_)
  }

  valid_T_star_b <- T_star_b[is.finite(T_star_b)]
  n_valid_T_star_b <- length(valid_T_star_b)
  if (n_valid_T_star_b == 0) {
    stop("Could not compute any valid bootstrap T* statistics. Cannot determine K_alpha. Check model fits and variance estimates.")
  }
  n_nonfinite_T <- B - n_valid_T_star_b
  if (n_nonfinite_T > 0) {
    add_warning(paste(n_nonfinite_T, "out of", B, "bootstrap T* statistics were non-finite.",
                      "This might be due to fitting errors or near-zero variance estimates.",
                      ifelse(n_nonfinite_T > B * 0.5, " Results may be unreliable.", "")), "T*b Calculation")
  }
  K_alpha <- stats::quantile(valid_T_star_b, probs = (1 - alpha), type = 8, na.rm = TRUE)

  # --- Construct Confidence Intervals ---
  if (verbose) message("Constructing confidence intervals...")
  results_list <- list()
  idx <- 1
  psi_ci_threshold <- 1e-12
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
        if (psi_hat > psi_ci_threshold) {
          se_nqj <- sqrt(psi_hat / n)
          margin <- K_alpha * se_nqj
          lower <- theta_qj - margin
          upper <- theta_qj + margin
        } else {
          lower <- theta_qj
          upper <- theta_qj
          se_nqj <- 0
          if(is.finite(theta_qj)){
            add_warning(paste("Psi_hat for '", coeff_name, "' in model '", q_name,
                              "' is near zero (", format(psi_hat, digits=3, scientific=TRUE),
                              "). Interval set to point estimate.", sep = ""), "CI Construction")
          }
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
  if (nrow(intervals_df) > 0) {
    intervals_df <- intervals_df[order(intervals_df$model_id, intervals_df$coefficient_name), ]
  }
  rownames(intervals_df) <- NULL

  # --- Output ---
  output <- list(
    intervals = intervals_df,
    K_alpha = as.numeric(K_alpha),
    alpha = alpha,
    B = B,
    n_valid_T_star_b = n_valid_T_star_b,
    T_star_b = T_star_b,
    bootstrap_method = bootstrap_method,
    warnings_list = warnings_list,
    valid_bootstrap_counts = valid_bs_counts,
    n_bootstrap_errors = total_bootstrap_errors
  )
  class(output) <- c("simultaneous_ci_result", "list")

  if (verbose) message("Done.")
  return(output)
}

