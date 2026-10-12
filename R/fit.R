# Stan orchestration. Source through R/load.R so .hsgp_root points to this
# directory. No compilation or sampling happens merely by sourcing this file.
#
# Compile in a separate directory. FLINT supplies the interval arithmetic used
# to bound the threshold throughout the shape interval. Returns a CmdStanModel
# for reuse through sample_hsgp(stan_model=...). constraint selects the Stan
# program, not the sign, which is handled by preparation and reconstruction.
# CmdStanR discovers CmdStan or uses cmdstanr::set_cmdstan_path(). Header-only
# edits require force_recompile=TRUE because make may not notice them.
build_hsgp <- function(constraint = c("monotone", "convex"),
                       build_dir = file.path(.hsgp_root, ".build"),
                       force_recompile = FALSE) {
  constraint <- match.arg(constraint)
  if (!requireNamespace("cmdstanr", quietly = TRUE)) stop("Install cmdstanr first.")
  pkg_config <- function(option) {
    value <- system2("pkg-config", c(option, "flint"), stdout = TRUE)
    if (!is.null(attr(value, "status"))) stop("FLINT 3 and pkg-config are required.")
    paste(value, collapse = " ")
  }
  flags <- pkg_config("--cflags")
  libs <- pkg_config("--libs")
  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)
  source <- file.path(.hsgp_root, "stan", paste0("hsgp_", constraint, "_marginal.stan"))
  target <- file.path(build_dir, basename(source))
  if (!file.exists(target) || !identical(readLines(source), readLines(target)))
    file.copy(source, target, overwrite = TRUE)
  header <- if (constraint == "monotone") "hsgp_global_stan_functions.hpp" else
    "hsgp_convex_stan_functions.hpp"
  # allow-undefined links the declarations in Stan to the supplied C++ header.
  # Suppressing FMA contraction preserves the interval arithmetic assumptions.
  cmdstanr::cmdstan_model(target,
    user_header = file.path(.hsgp_root, "stan", header),
    stanc_options = list("allow-undefined" = TRUE),
    cpp_options = list(CXXFLAGS = paste("-ffp-contract=off", flags), LDLIBS = libs),
    force_recompile = force_recompile)
}

# Sample the marginal coordinates and reconstruct the separated coefficient.
# Coefficients are returned in the original response units and sign convention.
# model comes from prepare_hsgp(). initial is a length-m vector in whitened
# coordinates, with zero used by default. initial_metric is Stan's dense
# inverse metric in those same coordinates (a covariance-shaped matrix).
# fixed_initial_metric disables metric adaptation while retaining step-size
# adaptation. warmup and samples are iteration counts per chain.
#
# Returns model, a draws-by-(m+1) coefficient matrix, chain diagnostics, the
# Stan summary, and the CmdStanMCMC object. summary uses Stan's internal sign
# for coefficient means/quantiles. predict_hsgp() uses the restored sign.
# fit.rds saves the numerical result without the live CmdStanMCMC object.
sample_hsgp <- function(model, output_dir, stan_model = NULL,
                        chains = 4L, parallel_chains = chains,
                        warmup = 1000L, samples = 1000L, seed = 124001L,
                        initial = NULL, initial_metric = NULL,
                        fixed_initial_metric = FALSE, max_treedepth = 12L,
                        adapt_delta = .95, stepsize = .01, refresh = 500L) {
  if (!inherits(model, "hsgp_model")) stop("Prepare a model with prepare_hsgp().")
  if (is.null(stan_model)) stan_model <- build_hsgp(model$constraint)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  parameter <- if (model$constraint == "monotone") "raw_z" else "z"
  init <- if (is.null(initial)) 0 else
    rep(list(setNames(list(as.numeric(initial)), parameter)), chains)
  options <- list(data = model$stan_data, seed = seed, chains = chains,
    parallel_chains = parallel_chains, iter_warmup = warmup,
    iter_sampling = samples, adapt_delta = adapt_delta,
    max_treedepth = max_treedepth, metric = "dense_e", step_size = stepsize,
    init = init, refresh = refresh, sig_figs = 18,
    output_dir = normalizePath(output_dir))
  if (!is.null(initial_metric)) options$inv_metric <- initial_metric
  if (fixed_initial_metric) {
    if (is.null(initial_metric)) stop("A fixed initial metric must be supplied.")
    # All warm-up is the final step-size adaptation stage. With no covariance
    # adaptation window, the supplied metric remains fixed, as in the lidar fit.
    options$init_buffer <- 0L
    options$window <- 0L
    options$term_buffer <- warmup
  }
  stan_fit <- do.call(stan_model$sample, options)
  # Generated quantities reconstruct one gamma per retained marginal draw.
  # Undo the response sign exactly once and remove posterior's draw attributes
  # so downstream linear algebra uses a plain coefficient matrix.
  coefficients <- as.matrix(stan_fit$draws("coefficients", format = "matrix")) * model$sign
  class(coefficients) <- c("matrix", "array")
  colnames(coefficients) <- paste0("a", 0:model$basis$rank)
  fit <- list(model = model, coefficients = coefficients,
    diagnostics = stan_fit$diagnostic_summary(),
    summary = stan_fit$summary(c(parameter, "coefficients")), stan_fit = stan_fit)
  saveRDS(fit[c("model", "coefficients", "diagnostics", "summary")],
          file.path(output_dir, "fit.rds"))
  write.csv(fit$summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
  print(fit$diagnostics)
  fit
}

# Descriptive aliases use the constraint stored in model to choose Stan.
sample_monotone_hsgp <- sample_hsgp
sample_convex_hsgp <- sample_hsgp

# Independent coefficient draws before imposing a shape constraint. Return
# the same model/coefficients pair accepted by predict_hsgp(). The factor F
# satisfies V=F F', so row-vector draws multiply by t(F), not by F.
# The seed is set explicitly, as in the application scripts.
draw_ordinary_hsgp <- function(model, samples = 20000L, seed = 20264000L) {
  set.seed(seed)
  p <- length(model$gaussian$mean)
  coefficients <- sweep(matrix(rnorm(samples * p), samples, p) %*%
    t(model$gaussian$factor), 2L, model$gaussian$mean, "+")
  list(model = model, coefficients = coefficients)
}
