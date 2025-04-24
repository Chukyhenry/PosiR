#' Plot Simultaneous Confidence Intervals
#'
#' Flexible plot method for objects returned by `simultaneous_ci`. Creates a
#' dotchart showing point estimates and confidence intervals, with options for
#' subsetting, scaling, rotation, and customization.
#'
#' @details
#' This function uses base R graphics. Layout, especially margin calculations
#' for labels, relies on approximations (e.g., using standard character heights/widths)
#' and may vary depending on the graphics device, device size, font, and specific
#' characters used in labels (including non-ASCII characters). For complex layouts
#' or guaranteed precision, consider using packages based on the 'grid' system
#' (like ggplot2).
#'
#' Faceting is not supported by this method.
#'
#' @param x An object of class `simultaneous_ci_result`.
#' @param y Ignored.
#' @param subset_pars Optional vector of coefficient names to include in the plot.
#'                    If NULL (default), all coefficients are plotted.
#' @param log.scale Logical. If TRUE, attempts to plot the estimates and intervals
#'                  on a logarithmic scale. Intervals crossing zero or non-positive
#'                  estimates/bounds will be omitted with a warning. Defaults to FALSE.
#' @param cex Numeric character expansion factor for estimate points (must be positive).
#'            Passed to `points`. Defaults to 0.8.
#' @param cex.labels Numeric character expansion factor for y-axis labels (must be positive).
#'                   Affects margin calculation. Defaults to 0.8.
#' @param las.labels Numeric value for label orientation relative to the axis (must be 0, 1, 2, or 3).
#'                   Passed to `Axis(side = 2, ...)`. Defaults to 1 (horizontal).
#'                   Use 2 (perpendicular) for long labels that might overlap.
#' @param pch Plotting character for estimate points. Passed to `points`.
#' @param col.estimate Color for the estimate points.
#' @param col.ci Color for the confidence interval segments.
#' @param col.ref Color for the reference line(s).
#' @param ref.line.pos Position(s) for vertical reference line(s). Defaults to 0.
#'                     Set to NULL to disable. Non-positive values are ignored if `log.scale=TRUE`.
#' @param lty.ref Line type for the reference line(s).
#' @param main Plot title.
#' @param xlab Label for the x-axis (estimates/intervals). If `log.scale=TRUE`,
#'             "Log Estimate" is appended unless `xlab` is explicitly set to NA.
#' @param label.trim Optional integer. If specified, trims labels longer than this
#'                   value using `substr` and adds "...". Helps with very long labels.
#' @param ... Additional arguments (currently ignored by this implementation but kept for S3 consistency).
#'
#' @return Invisibly returns a list containing:
#'   \item{ycoords}{A named numeric vector mapping labels to their y-axis coordinates.}
#'   \item{xlim}{The calculated x-axis limits.}
#'   \item{ylim}{The calculated y-axis limits.}
#'   Returns `invisible(NULL)` if no plot is generated.
#'
#' @method plot simultaneous_ci_result
#' @export
#' @importFrom graphics plot.new plot.window segments points Axis box title abline par strwidth grconvertX grconvertY
#' @importFrom stats setNames
#'
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
