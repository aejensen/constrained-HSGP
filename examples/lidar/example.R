# Mercury lidar preparation and summaries. See data/README.md for the data
# source. No optimization or sampling occurs when this file is sourced.
#
# Recorded parameter estimates. Predictor distances are in meters, while
# model x=(range_m-390)/330. Response offset/range refer to the negative log
# ratio. noise_log_sd describes the STANDARDIZED response noise, whereas
# amplitude_physical is in original response units and lengthscale_m in meters.
lidar_settings <- function() {
  list(rank = 128L, half_width = 2, center = 0.5, nu = 3.5,
    range_offset = 390, range_scale = 330,
    response_offset = 0.291156133819005, response_scale = 0.97646067,
    amplitude_physical = 0.308590337979769,
    lengthscale_m = 95.3267504759122,
    noise_log_sd = c(-2.85634865071989, 2.33025210851442, -0.305386704679671),
    intercept_mean_standardized = 0, intercept_sd_standardized = 2,
    absorption_factor = 1.6e-5, threshold_absolute_tolerance = 1e-12,
    max_precision_bits = 2048L)
}

# Read every released measurement in distance order and reverse the log ratio.
# The observed data frame retains both response signs and physical distances.
read_lidar_observations <- function(example_dir) {
  data <- read.csv(file.path(example_dir, "data", "lidar.csv"))
  if (!identical(names(data), c("range_m", "logratio")) || nrow(data) != 221L ||
      any(!is.finite(as.matrix(data))) || any(diff(data$range_m) <= 0) ||
      !identical(range(data$range_m), c(390L, 720L)))
    stop("The bundled example requires all 221 released lidar observations in order.")
  data$negative_logratio <- -data$logratio
  data
}

# Prepare the increasing model in standardized response units. Return the
# model together with settings and observations needed to undo that scaling.
prepare_lidar_example <- function(example_dir) {
  settings <- lidar_settings()
  observed <- read_lidar_observations(example_dir)
  x <- (observed$range_m - settings$range_offset) / settings$range_scale
  y <- (observed$negative_logratio - settings$response_offset) / settings$response_scale
  centered <- x - settings$center
  noise_sd <- exp(settings$noise_log_sd[1L] +
    settings$noise_log_sd[2L] * centered + settings$noise_log_sd[3L] * centered^2)
  # Scaling belongs to the example. The public API works in the units supplied.
  model <- prepare_monotone_hsgp(x, y, noise_sd = noise_sd,
    kernel = "matern", nu = settings$nu,
    amplitude = settings$amplitude_physical / settings$response_scale,
    lengthscale = settings$lengthscale_m / settings$range_scale,
    rank = settings$rank, interval = c(0, 1), half_width = settings$half_width,
    center = settings$center, intercept_mean = settings$intercept_mean_standardized,
    intercept_sd = settings$intercept_sd_standardized, direction = "increasing",
    threshold_absolute_tolerance = settings$threshold_absolute_tolerance,
    max_precision_bits = settings$max_precision_bits)
  list(model = model, settings = settings, observed = observed)
}

# Pointwise latent means and 95% equal-tail intervals on a common plotting
# grid. Without fit, return only the analytic ordinary Gaussian summaries.
# With fit, append summaries from the constrained coefficient draws.
lidar_plot_summaries <- function(prepared, fit = NULL) {
  settings <- prepared$settings
  observed_x <- (prepared$observed$range_m - settings$range_offset) / settings$range_scale
  grid <- sort(unique(c(seq(0, 1, length.out = 1001L), observed_x,
    head(observed_x, -1L) + diff(observed_x) / 2)))
  physical_summaries <- function(level, derivative, method) {
    # Chain rule: d/d(range_m)=(1/330)*d/dx. Restore response scale, then
    # divide by the DIAL absorption factor to obtain concentration in ng/m^3.
    concentration_scale <- settings$response_scale /
      (settings$range_scale * settings$absorption_factor)
    data.frame(method = method, x = level$x,
      range_m = settings$range_offset + settings$range_scale * level$x,
      mean = settings$response_offset + settings$response_scale * level$mean,
      lower = settings$response_offset + settings$response_scale * level$lower,
      upper = settings$response_offset + settings$response_scale * level$upper,
      mercury_ng_m3_mean = concentration_scale * derivative$mean,
      mercury_ng_m3_lower = concentration_scale * derivative$lower,
      mercury_ng_m3_upper = concentration_scale * derivative$upper)
  }
  ordinary <- physical_summaries(
    ordinary_hsgp(prepared$model, grid, derivative = 0L),
    ordinary_hsgp(prepared$model, grid, derivative = 1L), "ordinary_hsgp")
  if (is.null(fit)) return(ordinary)
  constrained <- physical_summaries(
    predict_hsgp(fit, grid, derivative = 0L),
    predict_hsgp(fit, grid, derivative = 1L), "global_monotone_hsgp")
  rbind(ordinary, constrained)
}

# Compare five residual discrepancies with one Gaussian replicate per curve.
# Rows are posterior draws and columns are released distances. Computing
# replicated standardized residuals directly as N(0,1) is equivalent to
# generating y_rep=f+noise first. Tail fractions reuse the fitting data.
lidar_predictive_checks <- function(prepared, fit, seed = 20264100L) {
  model <- prepared$model
  latent <- fit$coefficients %*% t(hsgp_basis(model$observations$x,
    model$basis$rank, model$basis$half_width, model$basis$center))
  residual <- sweep(sweep(-latent, 2, model$observations$y, "+"),
                    2, model$observations$noise_sd, "/")
  set.seed(seed)
  replicated <- matrix(rnorm(length(residual)), nrow(residual))
  left <- prepared$observed$range_m <= 555
  # Adjacency is in recorded distance order, with varying physical spacing.
  discrepancies <- function(e) cbind(
    mean_square = rowMeans(e^2), maximum = apply(abs(e), 1, max),
    adjacent_correlation = apply(e, 1, function(v) abs(cor(head(v, -1), tail(v, -1)))),
    mean_left = rowMeans(e[, left, drop = FALSE]),
    mean_right = rowMeans(e[, !left, drop = FALSE]))
  observed <- discrepancies(residual)
  replicated <- discrepancies(replicated)
  data.frame(discrepancy = colnames(observed),
    lower_tail = colMeans(replicated <= observed), upper_tail = colMeans(replicated >= observed))
}

# Optional re-estimation of the covariance and noise parameters by the ordinary
# underlying-GP marginal likelihood, using the same four starts as the paper.
# Return amplitude and length scale in physical units, standardized-noise
# log coefficients, and optim's convergence code. This leaves the recorded
# settings unchanged, allowing comparison before choosing to replace them.
estimate_lidar_parameters <- function(observed) {
  x <- (observed$range_m - 390) / 330
  y <- observed$negative_logratio
  y <- (y - mean(y)) / diff(range(y))
  z <- x - .5
  noise_design <- cbind(1, z, z^2)
  distance <- abs(outer(x, x, "-"))
  objective <- function(par) {
    q <- sqrt(7) * distance / exp(par[2])
    covariance <- exp(2 * par[1] - q) * (1 + q + 2*q^2/5 + q^3/15)
    noise <- exp(drop(noise_design %*% par[3:5]))
    # The constant 4 is the intercept variance (standardized intercept SD 2).
    C <- tryCatch(chol(covariance + 4 + diag(noise^2)), error = function(e) NULL)
    if (is.null(C)) return(1e100)
    sum(log(diag(C))) + sum(forwardsolve(t(C), y)^2)/2 + length(y)*log(2*pi)/2
  }
  starts <- list(c(log(c(.5, .08, .08)), 0, 0), c(log(c(1, .25, .05)), 2, 0),
                 c(log(c(2, .6, .04)), 3, 2), c(log(c(.5, .2, .05)), 2, 2))
  fits <- lapply(starts, function(p) optim(p, objective, method = "L-BFGS-B",
    lower = c(log(c(.005, .003, 1e-5)), -8, -8),
    upper = c(log(c(30, 3, .5)), 8, 8), control = list(maxit = 250, factr = 1e7)))
  best <- fits[[which.min(vapply(fits, `[[`, 0., "value"))]]
  list(amplitude_physical = exp(best$par[1]) * diff(range(observed$negative_logratio)),
       lengthscale_m = exp(best$par[2]) * 330,
       noise_log_sd = best$par[3:5], convergence = best$convergence)
}
