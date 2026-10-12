# Leaf light-response preparation and summaries. Data attribution and column
# units are in data/README.md. Sourcing this file has no fitting side effects.
#
# Apply the source study's quality filters to leaf 1818, then remove the dark
# measurement (PARi<15). Return 13 rows in irradiance order, retaining Time
# for the predictive assessment of adjacent observations.
read_leaf_observations <- function(example_dir = "examples/leaf") {
  leaf <- read.csv(file.path(example_dir, "data", "leaf1818_raw.csv"))
  leaf <- leaf[is.na(leaf$QC), ]
  leaf <- leaf[leaf$Cond >= .05 & leaf$Ci >= 0 & leaf$Ci <= 2000, ]
  leaf <- leaf[abs(leaf$Tleaf - mean(leaf$Tleaf)) <= .9999 & leaf$PARi >= 15, ]
  leaf[order(leaf$PARi), ]
}

# Map incident irradiance q to x=(q-origin)/scale in [-1,1]. Assimilation
# stays in its recorded units. The returned origin/scale convert derivatives
# back to q. The listed estimates come from the rank-32 ordinary fit below,
# and are shared by the rank-24 ordinary and concave posterior calculations.
prepare_leaf_example <- function(example_dir = "examples/leaf") {
  observed <- read_leaf_observations(example_dir)
  origin <- mean(range(observed$PARi))
  scale <- diff(range(observed$PARi)) / 2
  x <- (observed$PARi - origin) / scale
  model <- prepare_convex_hsgp(x, observed$Photo, noise_sd = .370337659274625,
    kernel = "se", amplitude = 5.62517804591132, lengthscale = .303384075312839,
    rank = 24L, interval = c(-1, 1), center = 0, half_width = 2,
    intercept_mean = 0, intercept_sd = 20, direction = "concave",
    threshold_absolute_tolerance = 1e-8)
  list(model = model, observed = observed, origin = origin, scale = scale)
}

# A polynomial calculation supplies a starting point near the marginal mode.
# Stan evaluates its conservative continuum threshold during sampling.
# Return a length-m initial vector in the prepared model's whitened coordinates.
# polyroot() and Horner evaluation are used only for initialization here.
leaf_initial <- function(model) {
  d <- model$stan_data
  m <- d$m
  # Column j contains the power coefficients of U_{j-1}, in ascending order.
  U <- matrix(0, m, m)
  U[1, 1] <- 1
  U[2, 2] <- 2
  if (m > 2L) for (j in 3:m)
    U[, j] <- c(0, 2 * U[-m, j - 1L]) - U[, j - 2L]
  interval <- sort(cos(pi * (c(d$left, d$right) + d$L) / (2 * d$L)))
  polynomial <- function(a, w) {
    value <- numeric(length(w))
    for (coefficient in rev(a)) value <- value * w + coefficient
    value
  }
  objective <- function(z) {
    g <- drop(d$g_mean + d$g_factor %*% z)
    a <- drop(U[, 2:m, drop = FALSE] %*% ((2:m)^2 * g[-1]))
    # The curvature threshold is the maximum of Q_g at stationary roots
    # and both cosine endpoints. The intercept g[1] does not enter Q_g.
    roots <- polyroot(seq_len(m - 1L) * a[-1])
    roots <- Re(roots[abs(Im(roots)) < 1e-7 * (1 + abs(Re(roots)))])
    candidates <- c(interval, roots[roots > interval[1] & roots < interval[2]])
    threshold <- max(polynomial(a, candidates))
    mu <- d$b_mean + sum(d$b_loading * z)
    sum(z^2)/2 - pnorm((threshold - mu)/d$conditional_sd,
                       lower.tail = FALSE, log.p = TRUE)
  }
  optim(numeric(m), objective, method = "BFGS",
        control = list(maxit = 2000L, reltol = 1e-10))$par
}

# Optional ordinary marginal-likelihood fit at rank 32, before the rank-24 fit.
# Optimize amplitude, normalized length scale, and noise SD on the log scale.
# The independent intercept variance is 20^2. Return the three estimates,
# without changing the parameters used by prepare_leaf_example().
estimate_leaf_parameters <- function(observed = read_leaf_observations()) {
  x <- (observed$PARi - mean(range(observed$PARi))) / (diff(range(observed$PARi))/2)
  H <- hsgp_basis(x, 32L, 2, 0)
  objective <- function(par) {
    p <- exp(par)
    sd <- c(20, hsgp_spectral_sd("se", p[1], p[2], 32L, 2))
    K <- tcrossprod(sweep(H, 2, sd, "*")) + diag(p[3]^2, nrow(observed))
    C <- chol(K)
    sum(log(diag(C))) + sum(forwardsolve(t(C), observed$Photo)^2)/2
  }
  fit <- optim(log(c(15, .65, .6)), objective, method = "L-BFGS-B",
    lower = log(c(1, .15, .1)), upper = log(c(50, 3, 3)))
  setNames(exp(fit$par), c("amplitude", "lengthscale", "noise_sd"))
}

# Compact plot data from concave and ordinary coefficient draws. Derivatives
# are converted from x to irradiance q by scale^(-r). Only the curvature
# display receives the additional factor 10^6. Intervals are pointwise 95%.
leaf_plot_summaries <- function(prepared, fit, ordinary) {
  grid <- seq(-1, 1, length.out = 301L)
  quantities <- c("Net assimilation", "Marginal assimilation gain", "Curvature (scaled by 10^6)")
  rows <- list()
  for (r in 0:2) for (name in c("Concave", "Unconstrained")) {
    object <- if (name == "Concave") fit else ordinary
    s <- predict_hsgp(object, grid, derivative = r)
    multiplier <- if (r == 2) 1e6/prepared$scale^2 else prepared$scale^(-r)
    s[c("mean", "lower", "upper")] <- s[c("mean", "lower", "upper")] * multiplier
    s$sd <- NULL
    s$q <- prepared$origin + prepared$scale * grid
    s$model <- name
    s$quantity <- quantities[r + 1L]
    rows[[length(rows) + 1L]] <- s
  }
  do.call(rbind, rows)
}

# Reproduce the table's response and slope summaries at q=100,500,1500.
# Compute differences draw by draw before taking means and quantiles, thereby
# retaining posterior dependence between the two evaluation locations.
leaf_scalar_summaries <- function(prepared, fit, ordinary) {
  x <- (c(100, 500, 1500) - prepared$origin) / prepared$scale
  H <- hsgp_basis(x, 24L, 2, 0)
  D <- hsgp_basis(x, 24L, 2, 0, derivative = 1L) / prepared$scale
  rows <- lapply(list(Concave = fit, Unconstrained = ordinary), function(object) {
    A <- object$coefficients %*% t(H)
    slope <- object$coefficients %*% t(D)
    values <- cbind(A100 = A[,1], A500 = A[,2], A1500 = A[,3],
      slope100 = slope[,1], slope1500 = slope[,3], slope_change = slope[,1] - slope[,3],
      gain500_1500 = A[,3] - A[,2])
    data.frame(quantity = colnames(values), mean = colMeans(values),
      lower = apply(values, 2, quantile, .025), upper = apply(values, 2, quantile, .975))
  })
  do.call(rbind, Map(function(d, name) cbind(model = name, d), rows, names(rows)))
}

# For each posterior curve, compare observed standardized residuals with an
# independent N(0,1) replicate vector. The maximum and the absolute sum of
# adjacent residual products are evaluated on the same curve. Sequence order
# uses measurement Time, not the irradiance ordering used in plots.
# Return upper-tail fractions for each model and discrepancy.
leaf_predictive_checks <- function(prepared, fit, ordinary, seed = 20264100L) {
  set.seed(seed)
  model <- prepared$model
  H <- hsgp_basis(model$observations$x, model$basis$rank, model$basis$half_width, model$basis$center)
  sequence <- order(prepared$observed$Time)
  discrepancies <- function(e) cbind(maximum = apply(abs(e), 1, max),
    sequence = abs(rowSums(e[, head(sequence, -1), drop = FALSE] *
                          e[, tail(sequence, -1), drop = FALSE])))
  rows <- lapply(list(Concave = fit, Unconstrained = ordinary), function(object) {
    latent <- object$coefficients %*% t(H)
    residual <- sweep(sweep(-latent, 2, model$observations$y, "+"), 2, model$observations$noise_sd, "/")
    replicated <- matrix(rnorm(length(residual)), nrow(residual))
    colMeans(discrepancies(replicated) >= discrepancies(residual))
  })
  data.frame(model = names(rows), do.call(rbind, rows), row.names = NULL)
}
