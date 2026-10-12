// Convex HSGP posterior, integrating gamma=-a_1 out of the Gaussian law.
// R/model.R supplies that law after conditioning on observations. The free
// vector is g=(a0,a2,...,am), and g=g_mean+g_factor*z whitens its marginal.
// Concave fits use sign-reversed data in R and reverse coefficients on return.
functions {
  // External C++ functions in hsgp_convex_stan_functions.hpp. Their threshold
  // covers the full interval through the Chebyshev polynomial for curvature.
  real convex_threshold(vector g, data real left, data real right,
                        data real L, data real tolerance);
  real convex_log_survival(real alpha);
  real convex_log_quantile(real log_probability);
}
data {
  // The rank limit is a bound of the numerical polynomial implementation.
  int<lower=2,upper=32> m;
  vector[m] g_mean;
  matrix[m,m] g_factor;
  // b denotes gamma=-a1. Given z, it has mean b_mean+b_loading'z and
  // variance conditional_sd^2. The square g_factor need not be triangular.
  real b_mean;
  vector[m] b_loading;
  real<lower=0> conditional_sd;
  // R has subtracted the basis center from both interval endpoints.
  // L is the padded half-width. The sine factor is positive on (-L,L).
  real left;
  real right;
  real<lower=0> L;
  // Absolute width allowed for the threshold bounds, in coefficient units.
  real<lower=0> tolerance;
}
transformed data {
  if (!(left < right && left > -L && right < L)) reject("Invalid shape interval.");
  if (!(conditional_sd > 0 && tolerance > 0)) reject("Positive SD/tolerance required.");
}
parameters {
  // The shape event reweights this ordinary standard normal marginal.
  vector[m] z;
}
transformed parameters {
  vector[m] g = g_mean + g_factor * z;
  real mu = b_mean + dot_product(b_loading,z);
  // The threshold is max_w sum_{j=2}^m j^2*g_j*U_{j-1}(w), with the
  // first entry of g holding the intercept and its second entry holding a2.
  real tau_upper = convex_threshold(g,left,right,L,tolerance);
  real alpha = (tau_upper - mu) / conditional_sd;
  real log_survival = convex_log_survival(alpha);
}
model {
  z ~ std_normal();
  // Probability that b exceeds the threshold conditional on z.
  target += log_survival;
}
generated quantities {
  vector[m+1] coefficients;
  // eta is a standard normal draw truncated below alpha. The log-tail
  // inverse remains usable when the ordinary tail probability underflows.
  real eta = -convex_log_quantile(log(uniform_rng(0,1)) + log_survival);
  // This form avoids cancellation between a large mean and a tail draw.
  real b = tau_upper + conditional_sd * (eta - alpha);
  if (is_nan(eta) || is_inf(eta) || eta < alpha) reject("Invalid reconstruction.");
  // Inverse map to (a0,a1,...,am), still in the internal convex sign.
  coefficients[1] = g[1];
  coefficients[2] = -b;
  coefficients[3:(m+1)] = g[2:m];
}
