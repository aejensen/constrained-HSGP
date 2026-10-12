# Gaussian preparation for monotone/convex HSGPs and latent-curve summaries.
# Sourcing this file does not compile a model, sample, or change the RNG state.

# Input helpers retain only finite numerical values and enforce the declared
# scalar, integer, or vector shape. They do not rescale the observations.
.hsgp_scalar <- function(value, name, positive = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      (positive && value <= 0))
    stop(name, " must be a finite ", if (positive) "positive " else "",
         "number.", call. = FALSE)
  as.numeric(value)
}

.hsgp_integer <- function(value, name, lower, upper = .Machine$integer.max) {
  value <- .hsgp_scalar(value, name)
  if (value < lower || value > upper || value != floor(value))
    stop(name, " must be an integer from ", lower, " through ", upper, ".",
         call. = FALSE)
  as.integer(value)
}

.hsgp_vector <- function(value, name) {
  if (!is.numeric(value) || !is.null(dim(value)) || any(!is.finite(value)))
    stop(name, " must be a finite numeric vector.", call. = FALSE)
  as.numeric(value)
}

# Return ORIGINAL spectral coefficient SDs sqrt(S(j*pi/(2*half_width))).
# Supported kernels are "matern" and "se". Amplitude is the stationary
# process SD, and lengthscale and half_width use the same units as x.
# The Fourier convention is k(r) = integral exp(i*omega*r) S(omega) d omega/(2*pi).
# Any nu > 0 defines a finite Matérn expansion. The paper's underlying-process
# convergence conditions require nu > r for a restriction of derivative order r.
hsgp_spectral_sd <- function(kernel = "matern", amplitude, lengthscale, rank,
                             half_width, nu = 3.5) {
  kernel <- match.arg(kernel, c("matern", "se"))
  amplitude <- .hsgp_scalar(amplitude, "amplitude", positive = TRUE)
  lengthscale <- .hsgp_scalar(lengthscale, "lengthscale", positive = TRUE)
  rank <- .hsgp_integer(rank, "rank", 1L)
  half_width <- .hsgp_scalar(half_width, "half_width", positive = TRUE)
  omega <- seq_len(rank) * pi / (2 * half_width)
  if (any(!is.finite(omega)) || any(omega <= 0))
    stop("Basis frequencies are outside floating-point range.", call. = FALSE)
  if (kernel == "se") {
    log_spectrum <- 2 * log(amplitude) + log(lengthscale) + log(2 * pi) / 2 -
      (lengthscale * omega)^2 / 2
  } else {
    nu <- .hsgp_scalar(nu, "nu", positive = TRUE)
    log_kappa <- (log(2) + log(nu)) / 2 - log(lengthscale)
    log_omega <- log(omega)
    # Evaluate log(kappa^2+omega^2) without squaring large/small numbers.
    log_denominator <- 2 * pmax(log_kappa, log_omega) +
      log1p(exp(-2 * abs(log_kappa - log_omega)))
    log_spectrum <- 2 * log(amplitude) + log(2) + log(pi) / 2 +
      lgamma(nu + 0.5) - lgamma(nu) + 2 * nu * log_kappa -
      (nu + 0.5) * log_denominator
  }
  result <- exp(log_spectrum / 2)
  if (any(!is.finite(result)) || any(result <= 0))
    stop("Spectral SDs are outside floating-point range. Adjust the rank or scales.",
         call. = FALSE)
  result
}

# Return a location-by-coefficient Laplacian basis matrix. With intercept=TRUE,
# columns are (a0,a1,...,am), and the intercept column is zero for derivatives.
# phi_j(x)=sin(omega_j*(x-center+half_width))/sqrt(half_width), with
# omega_j=j*pi/(2*half_width). Derivatives multiply by omega_j^r and cycle
# through sine/cosine and their negatives. Derivative orders 0 through 3
# are available. Evaluation does not impose a
# shape constraint, whose guarantee applies only on the prepared interval.
hsgp_basis <- function(x, rank, half_width, center, derivative = 0L,
                       intercept = TRUE) {
  x <- .hsgp_vector(x, "x")
  rank <- .hsgp_integer(rank, "rank", 1L)
  half_width <- .hsgp_scalar(half_width, "half_width", positive = TRUE)
  center <- .hsgp_scalar(center, "center")
  derivative <- .hsgp_integer(derivative, "derivative", 0L, 3L)
  if (!is.logical(intercept) || length(intercept) != 1L || is.na(intercept))
    stop("intercept must be TRUE or FALSE.", call. = FALSE)
  omega <- seq_len(rank) * pi / (2 * half_width)
  angle <- outer(x - center + half_width, omega)
  if (any(!is.finite(omega)) || any(omega <= 0) || any(!is.finite(angle)))
    stop("Basis evaluation is outside floating-point range.", call. = FALSE)
  spectral <- switch(as.character(derivative),
    "0" = sin(angle), "1" = cos(angle), "2" = -sin(angle), "3" = -cos(angle))
  spectral <- sweep(spectral, 2L, omega^derivative / sqrt(half_width), "*")
  if (any(!is.finite(spectral)))
    stop("Basis derivatives are outside floating-point range.", call. = FALSE)
  colnames(spectral) <- paste0("a", seq_len(rank))
  if (intercept)
    cbind(a0 = rep(as.numeric(derivative == 0L), length(x)), spectral)
  else spectral
}

# Prepare a Gaussian regression posterior for an interval-wide sign constraint.
# No automatic rescaling occurs. x and interval use the original predictor units.
# Supply exactly one of noise_sd (a positive scalar or one SD per observation)
# and noise_cov (a symmetric positive-definite observation covariance matrix).
# With x=y=numeric(0), give interval explicitly to prepare the prior. No noise
# argument is then needed. The intercept always has strictly positive variance.
# A supplied positive spectral_sd vector overrides the built-in kernel and must
# have length rank. lengthscale is required only for a built-in kernel.
#
# direction chooses f'>=0, f'<=0, f''>=0, or f''<=0 through "increasing",
# "decreasing", "convex", or "concave". Each choice imposes only that sign.
# Monotonicity requires m>=4 and I inside center +/- L/2. Convexity requires
# m>=2 and I strictly inside center +/- L. The supplied solvers cap the ranks
# at 128 and 32 respectively, independently of the theoretical conditions.
# threshold_absolute_tolerance is an absolute threshold width in coefficient
# units. max_precision_bits controls the monotone solver's Arb precision.
# The convex solver uses its own 1024-bit cap in include/convex_threshold.hpp.
#
# The returned hsgp_model contains:
#   stan_data: inputs for the increasing or convex internal Stan model.
#   gaussian: mean and square factor of the ordinary posterior in ORIGINAL
#     response-sign coordinates, ordered (a0,a1,...,am). Its covariance is
#     tcrossprod(factor). The square factor need not be triangular.
#   basis: rank, center, half_width, and the interval carrying the shape constraint.
#   direction/sign: decreasing and concave fits use response -y and sign=-1.
#   observations/specification: original observations and fixed model parameters.
# Fitting wrappers must multiply reconstructed internal coefficients by sign
# exactly once before exposing them to predict_hsgp(). No parameter uncertainty
# is integrated by this preparation step.
prepare_hsgp <- function(x, y, noise_sd = NULL, noise_cov = NULL,
    kernel = "matern", amplitude = 1, lengthscale, rank = 32L,
    interval = range(x), half_width = diff(interval) * 2,
    center = mean(interval), intercept_mean = 0, intercept_sd = 1,
    direction = c("increasing", "decreasing", "convex", "concave"), spectral_sd = NULL,
    threshold_absolute_tolerance = 1e-10, max_precision_bits = 512L,
    nu = 3.5) {
  x <- .hsgp_vector(x, "x")
  y <- .hsgp_vector(y, "y")
  if (length(x) != length(y))
    stop("x and y must have equal length.", call. = FALSE)
  n <- length(x)
  direction <- match.arg(direction)
  constraint <- if (direction %in% c("increasing", "decreasing")) "monotone" else "convex"
  rank <- .hsgp_integer(rank, "rank",
    if (constraint == "monotone") 4L else 2L,
    if (constraint == "monotone") 128L else 32L)
  interval <- .hsgp_vector(interval, "interval")
  if (length(interval) != 2L || interval[1L] >= interval[2L])
    stop("interval must contain two strictly ordered endpoints.", call. = FALSE)
  half_width <- .hsgp_scalar(half_width, "half_width", positive = TRUE)
  center <- .hsgp_scalar(center, "center")
  if (constraint == "monotone" && max(abs(interval - center)) > half_width / 2)
    stop("interval must lie within center +/- half_width/2.", call. = FALSE)
  if (constraint == "convex" && max(abs(interval - center)) >= half_width)
    stop("The convexity interval must lie strictly inside center +/- half_width.", call. = FALSE)
  if (any(x < interval[1L] | x > interval[2L]))
    stop("Observation locations must lie in interval.", call. = FALSE)
  intercept_mean <- .hsgp_scalar(intercept_mean, "intercept_mean")
  intercept_sd <- .hsgp_scalar(intercept_sd, "intercept_sd", positive = TRUE)
  sign <- if (direction %in% c("increasing", "convex")) 1 else -1
  threshold_absolute_tolerance <- .hsgp_scalar(threshold_absolute_tolerance,
    "threshold_absolute_tolerance", positive = TRUE)
  max_precision_bits <- .hsgp_integer(max_precision_bits, "max_precision_bits",
    64L, 4096L)
  custom_spectrum <- !is.null(spectral_sd)
  if (custom_spectrum) {
    spectral_sd <- .hsgp_vector(spectral_sd, "spectral_sd")
    if (length(spectral_sd) != rank || any(spectral_sd <= 0))
      stop("spectral_sd must contain rank strictly positive SDs.", call. = FALSE)
  } else {
    kernel <- match.arg(kernel, c("matern", "se"))
    spectral_sd <- hsgp_spectral_sd(kernel, amplitude, lengthscale, rank,
                                    half_width, nu)
  }
  if (!is.null(noise_sd) && !is.null(noise_cov))
    stop("Supply only one of noise_sd and noise_cov.", call. = FALSE)
  if (n && is.null(noise_sd) && is.null(noise_cov))
    stop("Supply noise_sd or noise_cov for the observations.", call. = FALSE)
  noise_factor <- NULL
  if (!is.null(noise_cov)) {
    if (!is.matrix(noise_cov) || !is.numeric(noise_cov) ||
        !identical(dim(noise_cov), c(n, n)) || any(!is.finite(noise_cov)))
      stop("noise_cov must be a finite n-by-n numeric matrix.", call. = FALSE)
    if (n) {
      symmetry_error <- max(abs(noise_cov - t(noise_cov)))
      if (symmetry_error > 100 * .Machine$double.eps * max(abs(noise_cov)))
        stop("noise_cov must be symmetric.", call. = FALSE)
      noise_cov <- (noise_cov + t(noise_cov)) / 2
      noise_factor <- tryCatch(t(chol(noise_cov)), error = function(e)
        stop("noise_cov must be positive definite.", call. = FALSE))
    }
  } else if (!is.null(noise_sd)) {
    noise_sd <- .hsgp_vector(noise_sd, "noise_sd")
    if (!length(noise_sd) || !length(noise_sd) %in% c(1L, n) || any(noise_sd <= 0))
      stop("noise_sd must be positive, with one value or one per observation.",
           call. = FALSE)
    if (n && length(noise_sd) == 1L) noise_sd <- rep(noise_sd, n)
  }

  # Work first with prior-standardized coefficients u: a=prior_mean+diag(scales)u.
  # Whitening the observation errors turns conditioning into ridge regression
  # with precision I+H'H. QR avoids explicitly inverting that precision.
  scales <- c(intercept_sd, spectral_sd)
  p <- rank + 1L
  internal_mean <- c(sign * intercept_mean, numeric(rank))
  design <- hsgp_basis(x, rank, half_width, center)
  if (!n) {
    posterior_mean <- internal_mean
    posterior_factor <- diag(scales, nrow = p)
    precision_factor <- diag(p)
    precision_pivot <- seq_len(p)
  } else {
    standardized_design <- sweep(design, 2L, scales, "*")
    standardized_response <- sign * (y - intercept_mean)
    if (is.null(noise_factor)) {
      standardized_design <- sweep(standardized_design, 1L, noise_sd, "/")
      standardized_response <- standardized_response / noise_sd
    } else {
      standardized_design <- forwardsolve(noise_factor, standardized_design)
      standardized_response <- drop(forwardsolve(noise_factor, standardized_response))
    }
    if (any(!is.finite(c(standardized_design, standardized_response))))
      stop("Gaussian standardization is outside floating-point range.", call. = FALSE)
    # The augmented design includes I, so the prior remains proper even if
    # observation locations repeat or the design has dependent columns.
    regression_qr <- qr(rbind(diag(p), standardized_design), LAPACK = TRUE)
    precision_factor <- qr.R(regression_qr)
    precision_pivot <- regression_qr$pivot
    response <- qr.qty(regression_qr, c(numeric(p), standardized_response))
    white_mean <- numeric(p)
    white_mean[precision_pivot] <- backsolve(precision_factor, response[seq_len(p)])
    white_factor <- matrix(0, p, p)
    white_factor[precision_pivot, ] <- backsolve(precision_factor, diag(p))
    posterior_mean <- internal_mean + scales * white_mean
    posterior_factor <- sweep(white_factor, 1L, scales, "*")
  }

  # Separate gamma=-a2 for monotonicity, or gamma=-a1 for convexity.
  separated <- if (constraint == "monotone") 3L else 2L
  free <- setdiff(seq_len(p), separated)
  chart <- diag(p)[free, , drop = FALSE]
  if (constraint == "monotone") chart[which(free == 5L), separated] <- -1 / 4
  chart <- rbind(chart, replace(numeric(p), separated, -1))
  chart_mean <- drop(chart %*% posterior_mean)
  chart_factor <- chart %*% posterior_factor
  free_factor <- chart_factor[seq_len(rank), , drop = FALSE]
  # Project the separated coordinate on the span of the remaining ones.
  # QR of the transposed factor gives g=g_mean+B*z, z~N(0,I), and the
  # conditional mean of gamma is chart_mean[p]+conditional_loading'z.
  free_qr <- qr(t(free_factor), LAPACK = TRUE)
  free_r <- qr.R(free_qr)
  if (any(!is.finite(free_r)) || any(diag(free_r) == 0))
    stop("The marginal Gaussian factor is numerically singular.", call. = FALSE)
  free_q <- qr.Q(free_qr)
  marginal_factor <- matrix(0, rank, rank)
  marginal_factor[free_qr$pivot, ] <- t(free_r)
  conditional_loading <- drop(chart_factor[p, , drop = FALSE] %*% free_q)
  # Holding g fixed and increasing gamma by one moves a along witness.
  # Consequently Var(gamma | g)=1/(witness' V^{-1} witness), where V is the
  # ordinary posterior covariance. This avoids subtracting nearly equal
  # variances in a Schur complement. Scaling the norm avoids overflow.
  witness <- numeric(p)
  witness[separated] <- -1
  if (constraint == "monotone") witness[5L] <- -1 / 4
  precision_witness <- drop(precision_factor %*%
    (witness / scales)[precision_pivot])
  unit <- max(abs(precision_witness))
  conditional_sd <- (1 / unit) / sqrt(sum((precision_witness / unit)^2))
  if (!is.finite(conditional_sd) || conditional_sd <= 0 ||
      any(!is.finite(c(chart_mean, marginal_factor, conditional_loading))))
    stop("The Gaussian coordinates are outside floating-point range.", call. = FALSE)
  stan_data <- list(m = rank, n_intercept = 1L,
    left = interval[1L], right = interval[2L], center = center,
    half_width = half_width, intercept_mean = sign * intercept_mean,
    g_mean = chart_mean[seq_len(rank)], g_factor = marginal_factor,
    t_mean = chart_mean[p], t_loading = conditional_loading,
    conditional_sd = conditional_sd,
    threshold_absolute_tolerance = threshold_absolute_tolerance,
    max_precision_bits = max_precision_bits)
  # Both Stan models consume this same conditional Gaussian construction,
  # with their respective coordinate names. The convex solver uses a centered
  # domain, so subtract center from its endpoints without rescaling x.
  if (constraint == "convex") {
    stan_data <- list(m = rank, g_mean = chart_mean[seq_len(rank)],
      g_factor = marginal_factor, b_mean = chart_mean[p],
      b_loading = conditional_loading, conditional_sd = conditional_sd,
      left = interval[1L] - center, right = interval[2L] - center,
      L = half_width, tolerance = threshold_absolute_tolerance)
  }
  coefficient_names <- paste0("a", 0:rank)
  # Restore original response sign for public summaries. Stan's inputs stay
  # in the increasing/convex convention until sampled coefficients return.
  gaussian <- list(mean = setNames(sign * posterior_mean, coefficient_names),
    factor = sign * posterior_factor)
  rownames(gaussian$factor) <- coefficient_names
  structure(list(stan_data = stan_data, gaussian = gaussian,
    basis = list(rank = rank, center = center, half_width = half_width,
                 interval = interval), direction = direction, sign = sign, constraint = constraint,
    observations = list(x = x, y = y, noise_sd = noise_sd, noise_cov = noise_cov),
    specification = list(kernel = if (custom_spectrum) "custom" else kernel,
      amplitude = if (custom_spectrum) NULL else amplitude,
      lengthscale = if (custom_spectrum) NULL else lengthscale,
      nu = if (!custom_spectrum && kernel == "matern") nu else NULL,
      spectral_sd = spectral_sd, intercept_mean = intercept_mean,
      intercept_sd = intercept_sd)), class = c(paste0(constraint, "_hsgp_model"), "hsgp_model"))
}

# Shared prediction preparation: validate interval probabilities and evaluate
# the requested derivative basis. probs supplies the two pointwise tail cuts.
.hsgp_prediction_inputs <- function(model, new_x, derivative, probs) {
  if (!inherits(model, "hsgp_model"))
    stop("model must come from prepare_hsgp().", call. = FALSE)
  new_x <- .hsgp_vector(new_x, "new_x")
  probs <- .hsgp_vector(probs, "probs")
  if (length(probs) != 2L || probs[1L] <= 0 || probs[2L] >= 1 ||
      probs[1L] >= probs[2L])
    stop("probs must contain two increasing probabilities strictly between zero and one.",
         call. = FALSE)
  basis <- model$basis
  list(x = new_x, probs = probs, design = hsgp_basis(new_x, basis$rank,
    basis$half_width, basis$center, derivative))
}

# Pointwise latent Gaussian summaries before imposing the constraint. The returned
# data frame has x, mean, sd, lower, and upper, with quantiles at probs. These are
# conditional on the fixed model parameters and do not include measurement noise.
ordinary_hsgp <- function(model, new_x, derivative = 0L, probs = c(.025, .975)) {
  inputs <- .hsgp_prediction_inputs(model, new_x, derivative, probs)
  mean <- drop(inputs$design %*% model$gaussian$mean)
  sd <- sqrt(rowSums((inputs$design %*% model$gaussian$factor)^2))
  data.frame(x = inputs$x, mean = mean, sd = sd,
    lower = mean + stats::qnorm(probs[1L]) * sd,
    upper = mean + stats::qnorm(probs[2L]) * sd)
}

# Pointwise latent summaries from retained coefficient draws. object must contain
# model and a finite numeric coefficients matrix with at least two draw rows and
# columns (a0,a1,...,am), already in the ORIGINAL response-sign coordinates.
# Quantiles use R's type-7 empirical rule. Output columns match ordinary_hsgp().
# Derivatives are with respect to the x supplied to preparation. No scaling,
# measurement noise, or hyperparameter uncertainty is added here.
predict_hsgp <- function(object, new_x, derivative = 0L, probs = c(.025, .975)) {
  if (!is.list(object) || is.null(object$model))
    stop("object must contain model and coefficients.", call. = FALSE)
  inputs <- .hsgp_prediction_inputs(object$model, new_x, derivative, probs)
  coefficients <- object$coefficients
  if (!is.matrix(coefficients) || !is.numeric(coefficients) ||
      ncol(coefficients) != object$model$basis$rank + 1L ||
      nrow(coefficients) < 2L || any(!is.finite(coefficients)))
    stop("coefficients must have at least two finite draw rows and rank+1 columns.",
         call. = FALSE)
  if (!length(inputs$x))
    return(data.frame(x = numeric(), mean = numeric(), sd = numeric(),
                      lower = numeric(), upper = numeric()))
  values <- coefficients %*% t(inputs$design)
  intervals <- vapply(seq_len(ncol(values)), function(j)
    stats::quantile(values[, j], probs = probs, names = FALSE, type = 7L), numeric(2L))
  data.frame(x = inputs$x, mean = colMeans(values), sd = apply(values, 2L, stats::sd),
    lower = intervals[1L, ], upper = intervals[2L, ])
}

# Convenience entry point restricting direction to increasing/decreasing.
# All remaining arguments are passed unchanged to prepare_hsgp().
prepare_monotone_hsgp <- function(..., direction = c("increasing", "decreasing")) {
  prepare_hsgp(..., direction = match.arg(direction))
}

# Convenience entry point restricting direction to convex/concave.
# Concavity alone permits either sign of the first derivative.
prepare_convex_hsgp <- function(..., direction = c("convex", "concave")) {
  prepare_hsgp(..., direction = match.arg(direction))
}
