# Expectation propagation (EP) for the lidar sampler's starting point and
# inverse metric. The input is a model from prepare_monotone_hsgp(). A grid of
# derivative inequalities is approximated by Gaussian site factors in the
# whitened ordinary posterior, then projected into Stan's marginal coordinates.
#
# EP supplies an initial point and metric, not the posterior target.
# The public sampler still evaluates the full interval constraint at every step.
# grid_size chooses derivative sites, damping controls site updates, and
# max_sweeps caps iterations. Return initial (length m), inverse_metric (m by m),
# site grid, and convergence history. Stan calls this covariance-shaped matrix
# its inverse metric. Convergence here concerns EP, not posterior MCMC.
prepare_derivative_ep <- function(model, grid_size = 129L, damping = 0.5,
                                  max_sweeps = 500L) {
  if (!is.numeric(grid_size) || length(grid_size) != 1L ||
      !is.finite(grid_size) || grid_size < 2L || grid_size != as.integer(grid_size) ||
      !is.numeric(damping) || length(damping) != 1L ||
      !is.finite(damping) || damping <= 0 || damping > 1 ||
      !is.numeric(max_sweeps) || length(max_sweeps) != 1L ||
      !is.finite(max_sweeps) || max_sweeps < 1L || max_sweeps != as.integer(max_sweeps))
    stop("Use an integer grid size and sweep limit, and damping in (0, 1].")
  data <- model$stan_data
  gaussian <- model$gaussian
  rank <- as.integer(data$m)
  dimension <- rank + 1L
  if (data$n_intercept != 1L || length(gaussian$mean) != dimension ||
      !identical(dim(gaussian$factor), c(dimension, dimension)))
    stop("EP geometry requires the random-intercept Gaussian preparation.")
  # The preparation stores original-sign coefficients. Stan uses increasing ones.
  direction_sign <- model$sign
  if (length(direction_sign) != 1L || !direction_sign %in% c(-1, 1))
    stop("The prepared model must specify its monotonicity direction.")
  mean_coefficient <- direction_sign * gaussian$mean
  factor_coefficient <- direction_sign * gaussian$factor
  grid <- seq(data$left, data$right, length.out = grid_size)
  frequency <- seq_len(rank) * pi / (2 * data$half_width)
  derivative <- cbind(0, sweep(
    cos(outer(grid - data$center + data$half_width, frequency)),
    2L, frequency / sqrt(data$half_width), "*"))
  offset <- drop(derivative %*% mean_coefficient)
  loading <- derivative %*% factor_coefficient
  scale <- sqrt(rowSums(loading^2))
  if (any(!is.finite(scale)) || any(scale <= 0))
    stop("An EP derivative variance is nonpositive or nonfinite.")
  # Normalize each derivative by its marginal SD to keep site scales comparable.
  site_loading <- sweep(loading, 1L, scale, "/")
  site_offset <- offset / scale

  # Each site approximates one truncated-normal factor in full Gaussian whitening.
  site_precision <- site_natural <- numeric(grid_size)
  mean_white <- numeric(dimension)
  covariance_white <- diag(dimension)
  history <- vector("list", max_sweeps)
  started <- proc.time()[["elapsed"]]
  converged <- FALSE
  for (sweep_id in seq_len(max_sweeps)) {
    previous_mean <- mean_white
    previous_covariance <- covariance_white
    maximum_site_change <- 0
    indices <- if (sweep_id %% 2L) seq_len(grid_size) else rev(seq_len(grid_size))
    for (i in indices) {
      vector <- site_loading[i, ]
      covariance_loading <- drop(covariance_white %*% vector)
      variance <- sum(vector * covariance_loading)
      mean <- sum(vector * mean_white)
      # Remove the current site's Gaussian factor to form the cavity law.
      # Its truncated-normal moments define the replacement site parameters.
      cavity_denominator <- 1 - site_precision[i] * variance
      if (!is.finite(cavity_denominator) || cavity_denominator <= 0)
        stop("An EP cavity precision is nonpositive.")
      cavity_variance <- variance / cavity_denominator
      cavity_mean <- (mean - site_natural[i] * variance) / cavity_denominator
      cavity_sd <- sqrt(cavity_variance)
      standardized_mean <- (cavity_mean + site_offset[i]) / cavity_sd
      mills <- exp(dnorm(standardized_mean, log = TRUE) -
        pnorm(standardized_mean, log.p = TRUE))
      variance_ratio <- 1 - mills * (mills + standardized_mean)
      if (!is.finite(variance_ratio) || variance_ratio <= 0 || variance_ratio > 1 + 1e-12)
        stop("Truncated-normal EP moments are outside numerical range.")
      tilted_mean <- cavity_mean + cavity_sd * mills
      tilted_variance <- cavity_variance * min(1, variance_ratio)
      next_precision <- max(0, 1 / tilted_variance - 1 / cavity_variance)
      next_natural <- tilted_mean / tilted_variance - cavity_mean / cavity_variance
      # Damping tempers each change in natural parameters before the rank-one update.
      delta_precision <- damping * (next_precision - site_precision[i])
      delta_natural <- damping * (next_natural - site_natural[i])
      maximum_site_change <- max(maximum_site_change,
        abs(delta_precision) / (1 + abs(site_precision[i])),
        abs(delta_natural) / (1 + abs(site_natural[i])))
      denominator <- 1 + delta_precision * variance
      if (!is.finite(denominator) || denominator <= 0)
        stop("An EP covariance update is nonpositive.")
      mean_white <- mean_white + covariance_loading *
        (delta_natural - delta_precision * mean) / denominator
      covariance_white <- covariance_white -
        (delta_precision / denominator) * tcrossprod(covariance_loading)
      site_precision[i] <- site_precision[i] + delta_precision
      site_natural[i] <- site_natural[i] + delta_natural
    }
    # Re-form the Gaussian approximation after each sweep to limit update drift.
    precision <- diag(dimension) + crossprod(site_loading,
      sweep(site_loading, 1L, site_precision, "*"))
    precision_factor <- chol(precision)
    covariance_white <- chol2inv(precision_factor)
    natural <- drop(crossprod(site_loading, site_natural))
    mean_white <- drop(backsolve(precision_factor,
      forwardsolve(t(precision_factor), natural)))
    mean_change <- max(abs(mean_white - previous_mean))
    covariance_change <- max(abs(covariance_white - previous_covariance))
    history[[sweep_id]] <- data.frame(sweep = sweep_id,
      maximum_site_change = maximum_site_change,
      maximum_mean_change = mean_change,
      maximum_covariance_change = covariance_change)
    # These tolerances stop construction of a stable initialization. They do
    # not control the continuum threshold or the constrained posterior error.
    if (maximum_site_change < 1e-6 && mean_change < 1e-7 && covariance_change < 1e-7) {
      converged <- TRUE
      break
    }
  }
  if (!converged)
    stop("EP geometry did not reach its stability criterion. Increase max_sweeps or omit EP.")

  # Remove a_2 and replace a_4 by g_4 = a_4 - a_2/4, then use the actual
  # marginal factor. It need not be triangular because Gaussian preparation uses QR.
  free_indices <- setdiff(seq_len(dimension), 3L)
  free_chart <- diag(dimension)[free_indices, , drop = FALSE]
  free_chart[match(5L, free_indices), 3L] <- -1 / 4
  map <- solve(data$g_factor, free_chart %*% factor_coefficient)
  whitening_error <- max(abs(tcrossprod(map) - diag(rank)))
  mean_error <- max(abs(drop(free_chart %*% mean_coefficient) - data$g_mean))
  initial <- drop(map %*% mean_white)
  inverse_metric <- map %*% covariance_white %*% t(map)
  inverse_metric <- (inverse_metric + t(inverse_metric)) / 2
  if (any(!is.finite(initial)) || any(!is.finite(inverse_metric)) ||
      whitening_error > 1e-8 || mean_error > 1e-8)
    stop("The EP approximation and Stan coefficient chart are inconsistent.")
  invisible(chol(inverse_metric))
  list(initial = initial, inverse_metric = inverse_metric,
    history = do.call(rbind, history[seq_len(sweep_id)]),
    grid = grid, sweeps = sweep_id, damping = damping,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    whitening_error = whitening_error, chart_mean_error = mean_error,
    purpose = "Initialization and inverse metric only. The sampler uses the full interval constraint.")
}
