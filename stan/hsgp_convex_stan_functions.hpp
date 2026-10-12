#ifndef CONVEX_HSGP_STAN_HPP
#define CONVEX_HSGP_STAN_HPP

// Bridge between Stan's Eigen/autodiff types and the convexity threshold.
// Geometry and tolerances are data. Only the free coefficients g carry
// derivatives. The trailing ostream pointer is required by Stan's external
// function calling convention and is not used here.

#include <stan/math.hpp>
#include "../include/convex_threshold.hpp"
#include <ostream>

// Return the conservative upper threshold. Reverse mode attaches the
// envelope vector at the numerical maximizing location. This approximates
// a derivative of the exact threshold, not of the rounded upper endpoint.
// T is double for value calls or stan::math::var for gradient calls.
template<typename T,int R,int C,int O,int MR,int MC>
inline stan::return_type_t<T> convex_threshold(
  const Eigen::Matrix<T,R,C,O,MR,MC>& g,const double& left,const double& right,
  const double& L,const double& tolerance,std::ostream*) {
  std::vector<double> v(g.size());
  for(int i=0;i<g.size();++i) v[i]=stan::math::value_of(g(i));
  auto r=convex_hsgp::oracle(v.size(),left,right,L)(v,tolerance);
  if constexpr(stan::is_var<T>::value) {
    std::vector<stan::math::var> operands(g.size());
    for(int i=0;i<g.size();++i) operands[i]=g(i);
    return stan::math::precomputed_gradients(r.upper,operands,r.gradient);
  } else return r.upper;
}

// log P(Z>=alpha)=log Phi(-alpha). Its derivative is the negative Mills
// ratio -phi(alpha)/Phi(-alpha), computed without forming the small tail.
template<typename T>
inline stan::return_type_t<T> convex_log_survival(const T& alpha,std::ostream*) {
  double a=stan::math::value_of(alpha);
  double value=stan::math::std_normal_lcdf(-a);
  double gradient=-std::exp(-.5*a*a-.5*std::log(2*std::acos(-1.0))-value);
  if(!std::isfinite(value) || !std::isfinite(gradient)) throw std::domain_error("Normal survival outside numerical range.");
  if constexpr(stan::is_var<T>::value)
    return stan::math::precomputed_gradients(value,std::vector<stan::math::var>{alpha},std::vector<double>{gradient});
  else return value;
}

// Return q satisfying log Phi(q)=lp for generated-quantities inversion.
// Newton refinement uses a stable log-CDF residual. Failure is reported
// instead of clipping a draw to the truncation boundary.
inline double convex_log_quantile(const double& lp,std::ostream*) {
  double q=stan::math::std_normal_log_qf(lp);
  for(int i=0;i<4;++i) {
    double v=stan::math::std_normal_lcdf(q);
    double d=std::exp(-.5*q*q-.5*std::log(2*std::acos(-1.0))-v);
    q-=(v-lp)/d;
  }
  if(!std::isfinite(q) || std::abs(stan::math::std_normal_lcdf(q)-lp)>1e-10*std::max(1.0,std::abs(lp)))
    throw std::domain_error("Conditional normal quantile failed refinement.");
  return q;
}
#endif
