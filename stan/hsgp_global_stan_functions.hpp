#ifndef HSGP_GLOBAL_STAN_FUNCTIONS_HPP
#define HSGP_GLOBAL_STAN_FUNCTIONS_HPP

// Adapter only: posterior sampling is performed by Stan's sampler.
// The returned value is a certified upper threshold for the stored g values.
// The supplied gradient approximates the exact threshold's envelope gradient.
// It is NOT a certified derivative of the literal rounded upper endpoint.
// g follows the m spectral-mode slots, with g[1]=0 (mode 2). Geometry and
// tolerances are data and receive no derivatives. The final ostream pointer
// on each external function is Stan's required logging argument.
#include <stan/math.hpp>
#include "../include/hsgp_global_threshold_fast.hpp"
#include <atomic>
#include <cmath>
#include <limits>
#include <ostream>
#include <type_traits>
#include <vector>

namespace hsgp_global_stan_detail {
// Per-process counts, exported as generated quantities. A polynomial fallback
// is an alternative evaluation route. A failure denotes an unresolved
// threshold or normal-tail evaluation. Saturating increments avoid overflow.
inline std::atomic<int>& failure_count() {
  static std::atomic<int> count{0};
  return count;
}
inline std::atomic<int>& fallback_count() {
  static std::atomic<int> count{0};
  return count;
}
inline void record_fallback() {
  int current = fallback_count().load();
  while (current < std::numeric_limits<int>::max()
         && !fallback_count().compare_exchange_weak(current, current + 1)) {}
}
inline void record_failure() {
  int current = failure_count().load();
  while (current < std::numeric_limits<int>::max()
         && !failure_count().compare_exchange_weak(current, current + 1)) {}
}
}  // namespace hsgp_global_stan_detail

// Read the accumulated fallback count without changing it.
inline int hsgp_monotone_oracle_fallbacks(std::ostream*) {
  return hsgp_global_stan_detail::fallback_count().load();
}

// Read the accumulated threshold/normal-tail failure count.
inline int hsgp_monotone_oracle_failures(std::ostream*) {
  return hsgp_global_stan_detail::failure_count().load();
}

// log P(Z>=alpha)=log Phi(-alpha), with derivative -phi(alpha)/Phi(-alpha).
// Stan 2.37's lcdf primal is stable in these tails, but its separate
// reverse-mode derivative has a less accurate Mills-ratio approximation.
// Evaluate the analytic derivative from that stable primal instead.
template <typename T>
inline stan::return_type_t<T> hsgp_normal_log_survival(
    const T& alpha, std::ostream*) {
  static_assert(std::is_same<T, double>::value || stan::is_var<T>::value,
                "This function supports double and reverse-mode scalar var.");
  try {
    const double a = stan::math::value_of(alpha);
    const double value = stan::math::std_normal_lcdf(-a);
    const double derivative = -std::exp(
        -0.5 * a * a - 0.5 * std::log(2.0 * std::acos(-1.0)) - value);
    if (!std::isfinite(value) || !std::isfinite(derivative))
      throw std::domain_error("Normal log survival is outside its numerical range.");
    if constexpr (stan::is_var<T>::value) {
      return stan::math::precomputed_gradients(
          value, std::vector<stan::math::var>{alpha},
          std::vector<double>{derivative});
    } else {
      return value;
    }
  } catch (...) {
    hsgp_global_stan_detail::record_failure();
    throw;
  }
}

// This scalar function is used only in generated quantities. Newton
// refinement solves log Phi(q)=log_probability in very small log tails.
// The caller negates q to invert the survival function.
// A numerical failure is reported rather than clamping a conditional draw.
inline double hsgp_normal_log_qf(const double& log_probability, std::ostream*) {
  try {
    double quantile = stan::math::std_normal_log_qf(log_probability);
    if (!std::isfinite(quantile))
      throw std::domain_error("Conditional normal log probability is at a machine endpoint.");
    for (int iteration = 0; iteration < 3; ++iteration) {
      const double value = stan::math::std_normal_lcdf(quantile);
      const double derivative = std::exp(
          -0.5 * quantile * quantile - 0.5 * std::log(2.0 * std::acos(-1.0)) - value);
      if (!std::isfinite(value) || !std::isfinite(derivative) || derivative <= 0)
        throw std::domain_error("Normal log quantile refinement is unresolved.");
      quantile -= (value - log_probability) / derivative;
      if (!std::isfinite(quantile))
        throw std::domain_error("Normal log quantile refinement is nonfinite.");
    }
    const double residual = stan::math::std_normal_lcdf(quantile) - log_probability;
    const double tolerance = 32.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(log_probability));
    if (!std::isfinite(residual) || std::abs(residual) > tolerance)
      throw std::domain_error("Normal log quantile refinement missed its tolerance.");
    return quantile;
  } catch (...) {
    hsgp_global_stan_detail::record_failure();
    throw;
  }
}

// Strip autodiff values for the numerical solver, then attach its envelope
// vector with precomputed_gradients. The fast route requests an active
// location as well as a value enclosure and permits polynomial fallback.
template <typename T, int R, int C, int Options, int MR, int MC>
inline stan::return_type_t<T> hsgp_monotone_threshold(
    const Eigen::Matrix<T, R, C, Options, MR, MC>& g,
    const double& left, const double& right, const double& half_width,
    const double& center, const double& absolute_tolerance,
    const int& max_precision_bits, std::ostream*) {
  static_assert(std::is_same<T, double>::value || stan::is_var<T>::value,
                "This function supports double and reverse-mode scalar var.");
  std::vector<double> values(g.size());
  for (Eigen::Index j = 0; j < g.size(); ++j)
    values[j] = stan::math::value_of(g(j));
  hsgp_global::ThresholdOptions options;
  options.absolute_tolerance = absolute_tolerance;
  options.relative_tolerance = 0.0;
  options.initial_precision_bits = 64;
  options.max_precision_bits = max_precision_bits;
  hsgp_global::ThresholdResult result;
  hsgp_global::FastThresholdStats stats;
  try {
    result = hsgp_global::global_monotonicity_threshold_fast(
        values, left, right, half_width, center, options, &stats, true, true);
  } catch (...) {
    if (stats.used_fallback) hsgp_global_stan_detail::record_fallback();
    hsgp_global_stan_detail::record_failure();
    throw;
  }
  if (stats.used_fallback) hsgp_global_stan_detail::record_fallback();
  if constexpr (stan::is_var<T>::value) {
    std::vector<stan::math::var> operands(g.size());
    for (Eigen::Index j = 0; j < g.size(); ++j)
      operands[j] = g(j);
    return stan::math::precomputed_gradients(result.upper, operands,
                                             result.gradient);
  } else {
    return result.upper;
  }
}

#endif
