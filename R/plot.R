#' Plot Simultaneous Confidence Intervals
#'
#' Visualizes confidence intervals returned by `simultaneous_ci()` using base R graphics.
#' Estimates are shown as points with corresponding CI segments, grouped and labeled by
#' model and coefficient name. Supports customization for log scale, character sizes,
#' label trimming, and reference lines.
#'
#' @param x An object of class `simultaneous_ci_result`, typically returned by `simultaneous_ci()`.
#' @param y Ignored.
#' @param subset_pars Optional character vector. Coefficient names to subset the plot. Default: all.
#' @param log.scale Logical. Plot on logarithmic scale. Intervals crossing 0 or with nonpositive bounds are excluded.
#' @param cex Point size for estimates. Default = 0.8.
#' @param cex.labels Label size for y-axis. Default = 0.8.
#' @param las.labels Orientation of y-axis labels (0, 1, 2, or 3). Default = 1.
#' @param pch Plot character for point estimates. Default = 16.
#' @param col.estimate Color of point estimates. Default = "blue".
#' @param col.ci Color of confidence interval lines. Default = "darkgray".
#' @param col.ref Color of reference line(s). Default = "red".
#' @param ref.line.pos Position(s) for vertical reference line(s). Default = 0. Set to NULL to omit.
#' @param lty.ref Line type for reference lines. Default = 2 (dashed).
#' @param main Plot title. Default = "Simultaneous Confidence Intervals".
#' @param xlab X-axis label. If NULL and `log.scale = TRUE`, label defaults to "Log Estimate".
#' @param label.trim Integer. Trims long coefficient labels to this width (adds "..."). Optional.
#' @param ... Additional arguments passed for future use (currently ignored).
#'
#' @return Invisibly returns a list:
#'
#' - `ycoords`: Named vector of y-axis positions for each label
#' - `xlim`: Range of x-axis limits used
#' - `ylim`: Range of y-axis limits used
#'
#' If no valid intervals are available for plotting, returns `invisible(NULL)`.
#'
#' @method plot simultaneous_ci_result
#' @export
#' @importFrom graphics plot.new plot.window segments points Axis box title abline par strwidth grconvertX grconvertY
#' @importFrom stats setNames
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100*2), 100, 2, dimnames = list(NULL, c("X1", "X2")))
#' y <- 1 + X[,1] - X[,2] + rnorm(100)
#' res <- simultaneous_ci(X, y, list(mod = 1:3), B = 100, add_intercept = TRUE)
#' plot(res)
plot.simultaneous_ci_result <- function(x, y = NULL, subset_pars = NULL,
                                        log.scale = FALSE,
                                        cex = 0.8, cex.labels = 0.8, las.labels = 1, pch = 16,
                                        col.estimate = "blue", col.ci = "darkgray", col.ref = "red",
                                        ref.line.pos = 0, lty.ref = 2,
                                        main = "Simultaneous Confidence Intervals",
                                        xlab = NULL,
                                        label.trim = NULL,
                                        ...) {

  if (!inherits(x, "simultaneous_ci_result")) {
    stop("Object must be of class 'simultaneous_ci_result'")
  }
  if (!is.numeric(cex) || length(cex) != 1 || cex <= 0) {
    warning("Invalid 'cex' value; using default 0.8.", call. = FALSE); cex <- 0.8
  }
  if (!is.numeric(cex.labels) || length(cex.labels) != 1 || cex.labels <= 0) {
    warning("Invalid 'cex.labels' value; using default 0.8.", call. = FALSE); cex.labels <- 0.8
  }
  if (!is.numeric(las.labels) || length(las.labels) != 1 || !(las.labels %in% 0:3)) {
    warning("Invalid 'las.labels' value; must be 0, 1, 2, or 3. Using default 1.", call. = FALSE); las.labels <- 1
  }

  if (is.null(x$intervals) || nrow(x$intervals) == 0) {
    warning("No interval data found in the object to plot.")
    return(invisible(NULL))
  }
  plot_data <- x$intervals
  if (!is.null(subset_pars)) {
    plot_data <- plot_data[plot_data$coefficient_name %in% subset_pars, ]
    if (nrow(plot_data) == 0) {
      warning("No data remaining after subsetting parameters.")
      return(invisible(NULL))
    }
  }
  plot_data$label <- ifelse(is.na(plot_data$model_id),
                            as.character(plot_data$coefficient_name),
                            paste(plot_data$model_id, plot_data$coefficient_name, sep = ": "))
  if (!is.null(label.trim) && is.numeric(label.trim) && label.trim > 3) {
    needs_trim <- nchar(plot_data$label) > label.trim
    plot_data$label[needs_trim] <- paste0(substr(plot_data$label[needs_trim], 1, label.trim - 3), "...")
  }
  plot_data <- plot_data[order(plot_data$label), ]
  plot_data_valid <- plot_data[is.finite(plot_data$estimate), ]

  valid_ref_lines <- ref.line.pos
  if (log.scale) {
    if (is.null(xlab)) xlab <- "Log Estimate"
    log_valid_idx <- (plot_data_valid$estimate > 0) &
      (is.finite(plot_data_valid$lower) & plot_data_valid$lower > 0) &
      (is.finite(plot_data_valid$upper) & plot_data_valid$upper > 0)
    n_removed <- sum(!log_valid_idx)
    if (n_removed > 0) {
      warning(paste(n_removed, "intervals removed due to non-positive estimate or bounds on log scale."), call. = FALSE)
      plot_data_valid <- plot_data_valid[log_valid_idx, ]
    }
    if (nrow(plot_data_valid) == 0) {
      warning("No valid data remaining after log scale transformation.")
      return(invisible(NULL))
    }
    plot_data_valid$plot_est <- log(plot_data_valid$estimate)
    plot_data_valid$plot_lower <- log(plot_data_valid$lower)
    plot_data_valid$plot_upper <- log(plot_data_valid$upper)
    if (!is.null(valid_ref_lines)) {
      valid_ref_lines <- valid_ref_lines[valid_ref_lines > 0 & is.finite(valid_ref_lines)]
      if(length(valid_ref_lines) > 0) { valid_ref_lines <- log(valid_ref_lines) } else { valid_ref_lines <- NULL }
    }
  } else {
    if (is.null(xlab)) xlab <- "Estimate"
    plot_data_valid$plot_est <- plot_data_valid$estimate
    plot_data_valid$plot_lower <- plot_data_valid$lower
    plot_data_valid$plot_upper <- plot_data_valid$upper
    if (!is.null(valid_ref_lines)) {
      valid_ref_lines <- valid_ref_lines[is.finite(valid_ref_lines)]
      if(length(valid_ref_lines) == 0) valid_ref_lines <- NULL
    }
  }
  plot_data_valid <- plot_data_valid[is.finite(plot_data_valid$plot_est) &
                                       is.finite(plot_data_valid$plot_lower) &
                                       is.finite(plot_data_valid$plot_upper), ]
  if (nrow(plot_data_valid) == 0) {
    warning("No valid (non-NA) interval data to plot after filtering and transformation.")
    return(invisible(NULL))
  }

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)

  label_widths_in <- graphics::strwidth(plot_data_valid$label, units = "inches", cex = cex.labels)
  label_heights_in <- graphics::strheight("M", units = "inches", cex = cex.labels)
  if (las.labels == 1) {
    required_margin_in <- max(label_widths_in) + 0.25
    margin_side <- 2
  } else if (las.labels == 2) {
    required_margin_in <- max(label_heights_in) * 1.5
    margin_side <- 2
  } else {
    required_margin_in <- max(label_widths_in) + 0.25
    margin_side <- 2
  }
  current_margins_in <- graphics::par("mai")
  current_plot_width_in <- graphics::par("pin")[1]
  if (required_margin_in > current_margins_in[margin_side]) {
    new_mai <- current_margins_in
    new_mai[margin_side] <- required_margin_in
    if (new_mai[2] + new_mai[4] >= current_plot_width_in * 0.95) {
      warning("Labels may be too long for the plot width on this device. ",
              "Consider subsetting, 'label.trim', reducing 'cex.labels'/'las.labels', or resizing plot device.", call. = FALSE)
    }
    graphics::par(mai = new_mai)
  }

  n_items <- nrow(plot_data_valid)
  ylim <- c(0.5, n_items + 0.5)
  xlim <- range(c(plot_data_valid$plot_lower, plot_data_valid$plot_upper), na.rm = TRUE)
  if(anyNA(xlim) || !all(is.finite(xlim))) {
    xlim <- range(plot_data_valid$plot_est, na.rm = TRUE)
  }
  xlim_range <- diff(xlim)
  if(is.finite(xlim_range) && xlim_range > 0) {
    xlim_buffer <- xlim_range * 0.05
  } else {
    xlim_buffer <- ifelse(is.finite(xlim[1]), abs(xlim[1] * 0.1), 1)
    if(xlim_buffer == 0) xlim_buffer <- 1
  }
  xlim <- xlim + c(-xlim_buffer, xlim_buffer)

  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = ylim, log = ifelse(log.scale, "x", ""))
  graphics::segments(plot_data_valid$plot_lower, 1:n_items,
                     plot_data_valid$plot_upper, 1:n_items,
                     col = col.ci)
  graphics::points(plot_data_valid$plot_est, 1:n_items, pch = pch, col = col.estimate, cex = cex)
  graphics::Axis(side = 1)
  graphics::Axis(side = 2, at = 1:n_items, labels = plot_data_valid$label, las = las.labels, cex.axis = cex.labels)
  graphics::box()
  graphics::title(main = main, xlab = xlab)
  if (!is.null(valid_ref_lines) && length(valid_ref_lines) > 0) {
    graphics::abline(v = valid_ref_lines, lty = lty.ref, col = col.ref)
  }

  ycoords <- stats::setNames(1:n_items, plot_data_valid$label)
  return_obj <- list(ycoords = ycoords, xlim = xlim, ylim = ylim)
  invisible(return_obj)
}
