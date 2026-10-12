// Monotone HSGP posterior in whitened marginal coordinates.
// R/model.R first conditions the Gaussian coefficient prior on the data and
// supplies the Gaussian law of (g,gamma), where gamma=-a_2 and
// g_4=a_4-a_2/4. The remaining g_j are the unchanged a_j, including a_0.
// Thus a_2=-gamma and a_4=g_4-gamma/4 on reconstruction.
// Decreasing fits enter through sign-reversed observations in R, and R
// reverses the reconstructed coefficients once on return.
//
// Original Gaussian posterior, with its increasing coordinate integrated out.
// The certified threshold covers the complete interval, with exact root fallback.
// U is a conservative numerical threshold. Its target is an inner event.
functions {
  // Definitions and reverse-mode derivatives are supplied by the user header
  // hsgp_global_stan_functions.hpp (compile with --allow-undefined).
  real hsgp_monotone_threshold(vector g, data real left, data real right,
                             data real half_width, data real center,
                             data real absolute_tolerance,
                             data int max_precision_bits);
  int hsgp_monotone_oracle_failures();
  int hsgp_monotone_oracle_fallbacks();
  real hsgp_normal_log_survival(real alpha);
  real hsgp_normal_log_qf(real log_probability);
}
data {
  // m is the number of sine modes. The R interface always sets n_intercept=1
  // for a random intercept, giving m remaining coordinates after integration.
  int<lower=4, upper=128> m;
  int<lower=0, upper=1> n_intercept;
  real left;
  real right;
  real center;
  real<lower=0> half_width;
  real intercept_mean;
  // g_coordinates = g_mean + g_factor*raw_z. g_factor is a square covariance
  // factor, not necessarily a Cholesky factor. The order is (a0,a1,a3,g4,...).
  vector[m - 1 + n_intercept] g_mean;
  matrix[m - 1 + n_intercept, m - 1 + n_intercept] g_factor;
  // gamma | raw_z is N(t_mean+t_loading'raw_z, conditional_sd^2).
  real t_mean;
  vector[m - 1 + n_intercept] t_loading;
  real<lower=0> conditional_sd;
  // Maximum allowed threshold enclosure width, in coefficient units.
  real<lower=0> threshold_absolute_tolerance;
  int<lower=64, upper=4096> max_precision_bits;
}
transformed data {
  if (left >= right || half_width <= 0
      || fmax(abs(left - center), abs(right - center)) > half_width / 2)
    reject("Require the whole shape interval within center +/- half_width / 2.");
  if (conditional_sd <= 0 || threshold_absolute_tolerance <= 0)
    reject("Conditional SD and threshold tolerance must be strictly positive.");
}
parameters {
  // Standard normal before weighting by the conditional shape probability.
  vector[m - 1 + n_intercept] raw_z;
}
transformed parameters {
  vector[m - 1 + n_intercept] g_coordinates = g_mean + g_factor * raw_z;
  vector[m] g;
  real threshold_upper;
  real conditional_t_mean = t_mean + dot_product(t_loading, raw_z);
  real alpha;
  real log_survival;
  {
    // The C++ threshold expects one slot per sine mode, with mode 2 set to
    // zero. The intercept has zero derivative and is omitted from this vector.
    int k = n_intercept + 1;
    for (j in 1:m) {
      if (j == 2) {
        g[j] = 0;
      } else {
        g[j] = g_coordinates[k];
        k += 1;
      }
    }
  }
  threshold_upper = hsgp_monotone_threshold(
    g, left, right, half_width, center,
    threshold_absolute_tolerance, max_precision_bits);
  alpha = (threshold_upper - conditional_t_mean) / conditional_sd;
  // Symmetry avoids the finite positive-tail cutoff in std_normal_lccdf.
  log_survival = hsgp_normal_log_survival(alpha);
}
model {
  raw_z ~ std_normal();
  // Integrating gamma over [threshold_upper,infinity) gives this survival
  // weight. There is no additional normalizing constant depending on raw_z.
  target += log_survival;
}
generated quantities {
  // Invert the conditional survival probability on the log scale. This
  // avoids subtracting an extreme normal CDF from one or underflowing its tail.
  real log_uniform = log(uniform_rng(0, 1));
  real eta = -hsgp_normal_log_qf(log_survival + log_uniform);
  real t;
  real intercept = intercept_mean;
  vector[m] a = g;
  vector[m + 1] coefficients;
  int threshold_oracle_failures;
  int threshold_oracle_fallbacks;
  if (is_nan(eta) || is_inf(eta) || eta < alpha)
    reject("Conditional Gaussian tail reconstruction failed.");
  // This form avoids cancelling the conditional mean against a large tail.
  t = threshold_upper + conditional_sd * (eta - alpha);
  a[2] = -t;
  a[4] = g[4] - t / 4;
  if (n_intercept == 1)
    intercept = g_coordinates[1];
  // Output order (a0,a1,...,am) matches hsgp_basis() in the R interface.
  coefficients[1] = intercept;
  coefficients[2:(m + 1)] = a;
  // Cumulative numerical counters for this Stan process. Falling back to
  // root isolation is normal. A failure denotes an unresolved evaluation.
  threshold_oracle_failures = hsgp_monotone_oracle_failures();
  threshold_oracle_fallbacks = hsgp_monotone_oracle_fallbacks();
}
