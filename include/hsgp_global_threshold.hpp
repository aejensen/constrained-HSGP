#ifndef HSGP_GLOBAL_THRESHOLD_HPP
#define HSGP_GLOBAL_THRESHOLD_HPP

// Exact-polynomial fallback for the monotonicity threshold.
// With theta=pi*(x-center+L)/(2L) and w=cos(theta), define
//   P(w)=sum_{j=1}^m j*g_j*T_j(w), D(w)=1+4*w^2-8*w^4.
// The returned bounds enclose tau_1(g)=max_w -P(w)/D(w). The physical
// interval lies in the central half of the padded domain, where D>0.
// Here g_j denotes a spectral coefficient, not the Gaussian intercept.
//
// Exhaustive continuum threshold for the fixed HSGP witness
// w_2=-1, w_4=-1/4. Stored doubles are interpreted as exact real numbers.
// The enclosure is validated. Its midpoint, active point, and envelope
// gradient are ordinary floating-point approximations. In particular,
// differentiating an enclosure endpoint does not give its exact derivative.
// This backend does not certify the law produced by a Stan sampler.

#include <flint/arb.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/arb_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <gmp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace hsgp_global {

// Stop once upper-lower <= absolute_tolerance + relative_tolerance*scale,
// where scale=max(abs(lower),abs(upper)). Tolerances have coefficient units.
// Precision is increased geometrically up to max_precision_bits.
struct ThresholdOptions {
  double absolute_tolerance = 1e-10;
  double relative_tolerance = 0.0;
  long initial_precision_bits = 64;
  long max_precision_bits = 512;
  std::size_t max_modes = 128;
};

// lower/upper enclose the exact maximum for the supplied binary doubles.
// estimate is their midpoint, active_y is a candidate in cosine coordinates,
// and gradient is the envelope vector there (g_2 has derivative zero).
// unique_active records strict separation from all competing candidates.
// The remaining fields describe root isolation and endpoint ambiguity.
struct ThresholdResult {
  double lower = 0.0;
  double upper = 0.0;
  double estimate = 0.0;
  double active_y = 0.0;
  std::vector<double> gradient;
  bool unique_active = false;
  long precision_bits = 0;
  long real_roots = 0;
  bool endpoint_overlap = false;
};

namespace detail {
// Small owning wrappers initialize and clear FLINT/Arb objects. Converting a
// double through mpq_set_d preserves its exact binary rational value.
struct Rational {
  fmpq_t value;
  Rational() { fmpq_init(value); }
  ~Rational() { fmpq_clear(value); }
  Rational(const Rational&) = delete;
  Rational& operator=(const Rational&) = delete;
  void set(double x) {
    mpq_t q;
    mpq_init(q);
    mpq_set_d(q, x);
    fmpq_set_mpq(value, q);
    mpq_clear(q);
  }
};
struct QPolynomial {
  fmpq_poly_t value;
  QPolynomial() { fmpq_poly_init(value); }
  ~QPolynomial() { fmpq_poly_clear(value); }
  QPolynomial(const QPolynomial&) = delete;
  QPolynomial& operator=(const QPolynomial&) = delete;
};
struct ZPolynomial {
  fmpz_poly_t value;
  ZPolynomial() { fmpz_poly_init(value); }
  ~ZPolynomial() { fmpz_poly_clear(value); }
  ZPolynomial(const ZPolynomial&) = delete;
  ZPolynomial& operator=(const ZPolynomial&) = delete;
};
struct Ball {
  arb_t value;
  Ball() { arb_init(value); }
  ~Ball() { arb_clear(value); }
  Ball(const Ball& other) { arb_init(value); arb_set(value, other.value); }
  Ball& operator=(const Ball& other) { arb_set(value, other.value); return *this; }
};
struct APolynomial {
  arb_poly_t value;
  APolynomial() { arb_poly_init(value); }
  ~APolynomial() { arb_poly_clear(value); }
};
struct Bounds {
  arf_t lower, upper;
  Bounds() { arf_init(lower); arf_init(upper); }
  ~Bounds() { arf_clear(lower); arf_clear(upper); }
  void set(const arb_t a, long bits) { arb_get_interval_arf(lower, upper, a, bits); }
};
struct Roots {
  arb_ptr values;
  slong size;
  explicit Roots(slong n) : values(_arb_vec_init(n)), size(n) {}
  ~Roots() { _arb_vec_clear(values, size); }
};
// A root overlapping a domain endpoint contributes an upper bound only.
// It can supply a lower bound only when membership is established.
struct Candidate {
  Ball y, objective;
  bool definitely_inside;
  Candidate(const arb_t point, bool inside) : definitely_inside(inside) {
    arb_set(y.value, point);
  }
};

// Return (endpoint-center+L)/(2L) exactly, before multiplying by pi.
inline void endpoint_ratio(fmpq_t result, double endpoint, double center,
                           double half_width) {
  Rational x, c, width, divisor;
  x.set(endpoint);
  c.set(center);
  width.set(half_width);
  fmpq_sub(result, x.value, c.value);
  fmpq_add(result, result, width.value);
  fmpq_mul_si(divisor.value, width.value, 2);
  fmpq_div(result, result, divisor.value);
}

// Construct P and D over rational coefficients. Stationary points satisfy
// P'D-PD'=0. Remove repeated factors before isolating the distinct real
// roots, retaining stationary points of every multiplicity.
inline void exact_polynomials(QPolynomial& p, QPolynomial& denominator,
                             ZPolynomial& stationary,
                             const std::vector<double>& g) {
  QPolynomial previous, current, next, term, dp, dd, product;
  Rational coefficient;
  fmpq_poly_set_coeff_si(previous.value, 0, 1);
  fmpq_poly_set_coeff_si(current.value, 1, 1);
  for (std::size_t j = 1; j <= g.size(); ++j) {
    coefficient.set(g[j - 1]);
    fmpq_mul_ui(coefficient.value, coefficient.value, static_cast<ulong>(j));
    fmpq_poly_scalar_mul_fmpq(term.value, current.value, coefficient.value);
    fmpq_poly_add(p.value, p.value, term.value);
    fmpq_poly_shift_left(next.value, current.value, 1);
    fmpq_poly_scalar_mul_si(next.value, next.value, 2);
    fmpq_poly_sub(next.value, next.value, previous.value);
    fmpq_poly_set(previous.value, current.value);
    fmpq_poly_set(current.value, next.value);
  }
  fmpq_poly_set_coeff_si(denominator.value, 0, 1);
  fmpq_poly_set_coeff_si(denominator.value, 2, 4);
  fmpq_poly_set_coeff_si(denominator.value, 4, -8);
  fmpq_poly_derivative(dp.value, p.value);
  fmpq_poly_derivative(dd.value, denominator.value);
  fmpq_poly_mul(dp.value, dp.value, denominator.value);
  fmpq_poly_mul(product.value, p.value, dd.value);
  fmpq_poly_sub(dp.value, dp.value, product.value);
  fmpq_poly_get_numerator(stationary.value, dp.value);
  if (fmpz_poly_degree(stationary.value) > 0) {
    ZPolynomial derivative, common, square_free;
    fmpz_poly_derivative(derivative.value, stationary.value);
    fmpz_poly_gcd(common.value, stationary.value, derivative.value);
    fmpz_poly_divexact(square_free.value, stationary.value, common.value);
    fmpz_poly_primitive_part(stationary.value, square_free.value);
  }
}

// Partial derivatives of -P(w)/D(w) with respect to the spectral g_j,
// evaluated at one candidate maximizer. At a tie this selects one active
// envelope vector and does not assert differentiability of the maximum.
inline void envelope_gradient(std::vector<double>& gradient, double y) {
  const double d = 1.0 + 4.0*y*y - 8.0*y*y*y*y;
  double previous = 1.0, current = y;
  for (std::size_t j = 1; j <= gradient.size(); ++j) {
    gradient[j - 1] = (j == 2) ? 0.0 : -static_cast<double>(j)*current/d;
    const double next = 2.0*y*current - previous;
    previous = current;
    current = next;
  }
}
}  // namespace detail

// g has one entry per original sine mode, with g[1] exactly zero.
// Physical coordinates have center c and padded half-width L. This routine
// requires m>=4 and L>=2*max(|left-c|,|right-c|). It fails closed if its
// finite precision/degree budget cannot resolve the requested enclosure.
inline ThresholdResult global_monotonicity_threshold(
    const std::vector<double>& g, double left, double right,
    double half_width, double center = 0.0,
    const ThresholdOptions& options = ThresholdOptions{}) {
  using namespace detail;
  if (g.size() < 4 || g.size() > options.max_modes || options.max_modes > 4096)
    throw std::invalid_argument("global threshold: unsupported mode count");
  for (double coefficient : g)
    if (!std::isfinite(coefficient))
      throw std::invalid_argument("global threshold: nonfinite coefficient");
  if (g[1] != 0.0)
    throw std::invalid_argument("global threshold: second free coefficient must be zero");
  if (!std::isfinite(left) || !std::isfinite(right) || !std::isfinite(center) ||
      !std::isfinite(half_width) || !(left < right) || !(half_width > 0.0))
    throw std::invalid_argument("global threshold: invalid physical interval");
  if (!std::isfinite(options.absolute_tolerance) ||
      !std::isfinite(options.relative_tolerance) ||
      options.absolute_tolerance < 0.0 || options.relative_tolerance < 0.0 ||
      !(options.absolute_tolerance > 0.0 || options.relative_tolerance > 0.0) ||
      options.initial_precision_bits < 32 ||
      options.max_precision_bits < options.initial_precision_bits ||
      options.max_precision_bits > 4096)
    throw std::invalid_argument("global threshold: invalid precision budget");
  Rational left_ratio, right_ratio;
  endpoint_ratio(left_ratio.value, left, center, half_width);
  endpoint_ratio(right_ratio.value, right, center, half_width);
  Rational quarter, three_quarters;
  fmpq_set_si(quarter.value, 1, 4);
  fmpq_set_si(three_quarters.value, 3, 4);
  if (fmpq_cmp(left_ratio.value, quarter.value) < 0 ||
      fmpq_cmp(right_ratio.value, three_quarters.value) > 0)
    throw std::invalid_argument("global threshold: padding must satisfy L>=2R");

  QPolynomial p, denominator;
  ZPolynomial stationary;
  exact_polynomials(p, denominator, stationary, g);
  const slong degree = std::max<slong>(0, fmpz_poly_degree(stationary.value));
  for (long bits = options.initial_precision_bits;;
       bits = std::min(options.max_precision_bits, bits * 2)) {
    // Cosine reverses endpoint order on theta in [pi/4,3*pi/4]. All real
    // stationary roots that might be in this interval join its endpoints.
    // Constant quotients need only the endpoint candidates.
    Ball y_min, y_max;
    arb_cos_pi_fmpq(y_min.value, right_ratio.value, bits);
    arb_cos_pi_fmpq(y_max.value, left_ratio.value, bits);
    Bounds domain_min, domain_max;
    domain_min.set(y_min.value, bits);
    domain_max.set(y_max.value, bits);
    std::vector<Candidate> candidates;
    candidates.emplace_back(y_min.value, true);
    candidates.emplace_back(y_max.value, true);
    Roots roots(degree);
    const slong number = degree > 0 ?
      arb_fmpz_poly_real_roots(roots.values, stationary.value, 0, bits) : 0;
    if (number < 0 || number > degree)
      throw std::runtime_error("global threshold: real-root isolation failed");
    bool overlap = false;
    for (slong i = 0; i < number; ++i) {
      Bounds root;
      root.set(roots.values + i, bits);
      if (arf_cmp(root.upper, domain_min.lower) < 0 ||
          arf_cmp(root.lower, domain_max.upper) > 0) continue;
      const bool inside = arf_cmp(root.lower, domain_min.upper) > 0 &&
                          arf_cmp(root.upper, domain_max.lower) < 0;
      overlap = overlap || !inside;
      candidates.emplace_back(roots.values + i, inside);
    }
    // Evaluate every candidate with ball arithmetic. A near-endpoint root
    // that might lie outside cannot increase the feasible lower bound.
    APolynomial ap, ad;
    arb_poly_set_fmpq_poly(ap.value, p.value, bits);
    arb_poly_set_fmpq_poly(ad.value, denominator.value, bits);
    bool finite = true;
    double lower = -std::numeric_limits<double>::infinity();
    double upper = lower;
    std::size_t active = 0;
    double highest_midpoint = lower;
    std::vector<double> lowers, uppers;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
      Ball divisor;
      arb_poly_evaluate(divisor.value, ad.value, candidates[i].y.value, bits);
      if (!arb_is_positive(divisor.value)) { finite = false; break; }
      arb_poly_evaluate(candidates[i].objective.value, ap.value,
                        candidates[i].y.value, bits);
      arb_neg(candidates[i].objective.value, candidates[i].objective.value);
      arb_div(candidates[i].objective.value, candidates[i].objective.value,
              divisor.value, bits);
      Bounds objective;
      objective.set(candidates[i].objective.value, bits);
      const double lo = arf_get_d(objective.lower, ARF_RND_FLOOR);
      const double hi = arf_get_d(objective.upper, ARF_RND_CEIL);
      if (!std::isfinite(lo) || !std::isfinite(hi)) { finite = false; break; }
      lowers.push_back(lo);
      uppers.push_back(hi);
      if (candidates[i].definitely_inside) lower = std::max(lower, lo);
      upper = std::max(upper, hi);
      const double mid = arf_get_d(arb_midref(candidates[i].objective.value),
                                   ARF_RND_NEAR);
      if (mid > highest_midpoint) { highest_midpoint = mid; active = i; }
    }
    // Enclosure width, rather than root-location width, controls stopping.
    // If it is unresolved, rebuild the endpoint and root balls at more bits.
    const long double width = static_cast<long double>(upper) - lower;
    const long double tolerance = options.absolute_tolerance +
      static_cast<long double>(options.relative_tolerance) *
      std::max(std::abs(lower), std::abs(upper));
    if (finite && width <= tolerance) {
      ThresholdResult result;
      result.lower = lower;
      result.upper = upper;
      result.estimate = lower/2.0 + upper/2.0;
      result.active_y = arf_get_d(arb_midref(candidates[active].y.value), ARF_RND_NEAR);
      result.gradient.resize(g.size());
      envelope_gradient(result.gradient, result.active_y);
      result.unique_active = candidates[active].definitely_inside &&
                             !fmpz_poly_is_zero(stationary.value);
      for (std::size_t i = 0; i < candidates.size(); ++i)
        if (i != active && !(lowers[active] > uppers[i])) result.unique_active = false;
      result.precision_bits = bits;
      result.real_roots = number;
      result.endpoint_overlap = overlap;
      return result;
    }
    if (bits == options.max_precision_bits) break;
  }
  throw std::runtime_error("global threshold: enclosure unresolved within precision budget");
}

}  // namespace hsgp_global
#endif
