#ifndef HSGP_GLOBAL_THRESHOLD_FAST_HPP
#define HSGP_GLOBAL_THRESHOLD_FAST_HPP

// Certified branch-and-bound on the entire continuum. The trigonometric
// objective is the same as in the exact-polynomial implementation. Arb
// encloses each point's sine and cosine. A stable recurrence-residual bound
// encloses the finite-mode sums, and Taylor bounds enclose every interval.
// The exact-polynomial implementation remains available as a fallback.
// For theta=pi*(x-center+L)/(2L), maximize
//   F(theta) = -sum_j j*g_j*cos(j*theta) / (-2*cos(2*theta)-cos(4*theta)).
// This is the rational Chebyshev quotient in hsgp_global_threshold.hpp.
// Fast calculations use theta coordinates. The returned active_y is cos(theta).
// IEEE round-to-nearest arithmetic is assumed, with no fast-math reassociation
// or fused multiply-add contraction (compile with -ffp-contract=off).
#include "hsgp_global_threshold.hpp"
#include <queue>

namespace hsgp_global {
namespace fast_detail {
// Expand each ordinary arithmetic result by one representable value in the
// required direction. Arb supplies rigorous starting trigonometric bounds.
inline double down(double x) {
  return std::nextafter(x, -std::numeric_limits<double>::infinity());
}
inline double up(double x) {
  return std::nextafter(x, std::numeric_limits<double>::infinity());
}
struct Interval {
  double lo, hi;
  Interval(double x = 0.0) : lo(x), hi(x) {}
  Interval(double a, double b) : lo(a), hi(b) {}
};
inline Interval operator+(Interval a, Interval b) {
  return {down(a.lo+b.lo), up(a.hi+b.hi)};
}
inline Interval operator-(Interval a, Interval b) {
  return {down(a.lo-b.hi), up(a.hi-b.lo)};
}
inline Interval operator-(Interval a) { return {-a.hi, -a.lo}; }
inline Interval operator*(Interval a, Interval b) {
  const double p[] = {a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
  return {down(*std::min_element(p,p+4)),up(*std::max_element(p,p+4))};
}
inline Interval operator/(Interval a, Interval b) {
  if (b.lo<=0.0 && b.hi>=0.0)
    throw std::runtime_error("fast threshold: denominator interval contains zero");
  return a * Interval(down(1.0/b.hi),up(1.0/b.lo));
}
inline double max_abs(Interval a) { return std::max(std::abs(a.lo),std::abs(a.hi)); }
inline Interval centered(double x, double radius) {
  return {down(x-radius),up(x+radius)};
}
// Convert an Arb ball to an outward-rounded double interval.
inline Interval ball_interval(const arb_t x, long bits) {
  detail::Bounds b;
  b.set(x,bits);
  return {arf_get_d(b.lower,ARF_RND_FLOOR),arf_get_d(b.upper,ARF_RND_CEIL)};
}
// Value, first derivative, and second derivative of F with respect to theta.
struct Point {
  double theta;
  Interval value, first, second;
};
// Cache coefficient-dependent derivative bounds for one threshold call.
// s[r] bounds sum_j j^(r+1)*abs(g_j), for r=0,...,3. curvature and third_bound
// bound |F''| and |F'''| on the domain. gamma is a floating-point
// summation-error factor, unrelated to the separated Gaussian coefficient.
struct Evaluator {
  const std::vector<double>& g;
  double s[4] = {0,0,0,0};
  double curvature, third_bound, gamma;
  long evaluations = 0,precise_evaluations=0,precise_bits;
  explicit Evaluator(const std::vector<double>& coefficients,long max_bits=128)
    : g(coefficients),precise_bits(std::min<long>(128,max_bits)) {
    for (std::size_t k=0;k<g.size();++k) {
      Interval term(std::abs(g[k]));
      for (int r=0;r<4;++r) {
        term=term*Interval(static_cast<double>(k+1));
        s[r]=(Interval(s[r])+term).hi;
      }
    }
    const double u=std::ldexp(1.0,-53);
    const Interval nu=Interval(static_cast<double>(g.size()+8))*Interval(u);
    gamma=(nu/(Interval(1.0)-nu)).hi;
    // On the central half of the padded domain D>=1. The bound .98 allows
    // the tiny enlargement of endpoints by interval arithmetic. Bounds for
    // |D'|, |D''|, |D'''| are 8, 24, and 80 respectively.
    const Interval d(.98),d2=d*d,d3=d2*d,d4=d3*d;
    curvature=(Interval(s[2])/d + Interval(16.0)*Interval(s[1])/d2
      +Interval(24.0)*Interval(s[0])/d2
      +Interval(128.0)*Interval(s[0])/d3).hi;
    third_bound=(Interval(s[3])/d+Interval(24.0)*Interval(s[2])/d2
      +Interval(72.0)*Interval(s[1])/d2
      +Interval(384.0)*Interval(s[1])/d3
      +Interval(80.0)*Interval(s[0])/d2
      +Interval(1152.0)*Interval(s[0])/d3
      +Interval(3072.0)*Interval(s[0])/d4).hi;
  }
  // Re-evaluate each harmonic directly in Arb when double recurrences cannot
  // locate a maximum or resolve its value to the requested tolerance.
  Point precise(double theta) {
    ++evaluations;++precise_evaluations;
    using detail::Ball;
    Ball input,angle,cosine,sine,coefficient,term,p,a,b,c2,c4,s2,s4;
    arb_set_d(input.value,theta);
    for(std::size_t k=0;k<g.size();++k) {
      const ulong j=static_cast<ulong>(k+1);
      arb_mul_ui(angle.value,input.value,j,precise_bits);
      arb_sin_cos(sine.value,cosine.value,angle.value,precise_bits);
      if(j==2) {arb_set(c2.value,cosine.value);arb_set(s2.value,sine.value);}
      if(j==4) {arb_set(c4.value,cosine.value);arb_set(s4.value,sine.value);}
      arb_set_d(coefficient.value,g[k]);
      arb_mul_ui(coefficient.value,coefficient.value,j,precise_bits);
      arb_mul(term.value,coefficient.value,cosine.value,precise_bits);
      arb_add(p.value,p.value,term.value,precise_bits);
      arb_mul_ui(coefficient.value,coefficient.value,j,precise_bits);
      arb_mul(term.value,coefficient.value,sine.value,precise_bits);
      arb_add(a.value,a.value,term.value,precise_bits);
      arb_mul_ui(coefficient.value,coefficient.value,j,precise_bits);
      arb_mul(term.value,coefficient.value,cosine.value,precise_bits);
      arb_add(b.value,b.value,term.value,precise_bits);
    }
    Ball d,d1,d2,dd,ddd,value,first,second,tmp;
    arb_mul_si(d.value,c2.value,-2,precise_bits);
    arb_sub(d.value,d.value,c4.value,precise_bits);
    if(!arb_is_positive(d.value))
      throw std::runtime_error("fast threshold: precise denominator unresolved");
    arb_add(d1.value,s2.value,s4.value,precise_bits);
    arb_mul_ui(d1.value,d1.value,4,precise_bits);
    arb_mul_ui(d2.value,c2.value,8,precise_bits);
    arb_mul_ui(tmp.value,c4.value,16,precise_bits);
    arb_add(d2.value,d2.value,tmp.value,precise_bits);
    arb_mul(dd.value,d.value,d.value,precise_bits);
    arb_mul(ddd.value,dd.value,d.value,precise_bits);
    arb_div(value.value,p.value,d.value,precise_bits);arb_neg(value.value,value.value);
    arb_div(first.value,a.value,d.value,precise_bits);
    arb_mul(tmp.value,p.value,d1.value,precise_bits);
    arb_div(tmp.value,tmp.value,dd.value,precise_bits);
    arb_add(first.value,first.value,tmp.value,precise_bits);
    arb_div(second.value,b.value,d.value,precise_bits);
    arb_mul(tmp.value,a.value,d1.value,precise_bits);
    arb_mul_si(tmp.value,tmp.value,-2,precise_bits);
    arb_div(tmp.value,tmp.value,dd.value,precise_bits);
    arb_add(second.value,second.value,tmp.value,precise_bits);
    arb_mul(tmp.value,p.value,d2.value,precise_bits);
    arb_div(tmp.value,tmp.value,dd.value,precise_bits);
    arb_add(second.value,second.value,tmp.value,precise_bits);
    arb_mul(tmp.value,d1.value,d1.value,precise_bits);
    arb_mul(tmp.value,tmp.value,p.value,precise_bits);
    arb_mul_si(tmp.value,tmp.value,-2,precise_bits);
    arb_div(tmp.value,tmp.value,ddd.value,precise_bits);
    arb_add(second.value,second.value,tmp.value,precise_bits);
    return {theta,ball_interval(value.value,precise_bits),
      ball_interval(first.value,precise_bits),ball_interval(second.value,precise_bits)};
  }
  // Fast point evaluation: T_j(cos(theta)) gives cos(j*theta), and
  // sin(theta)*U_{j-1}(cos(theta)) gives sin(j*theta). Recurrence and summation
  // errors are bounded explicitly before differentiating the quotient.
  Point operator()(double theta) {
    ++evaluations;
    detail::Ball input, cosine, sine;
    arb_set_d(input.value,theta);
    arb_sin_cos(sine.value,cosine.value,input.value,96);
    const Interval yc=ball_interval(cosine.value,96),sc=ball_interval(sine.value,96);
    const double y=arf_get_d(arb_midref(cosine.value),ARF_RND_NEAR);
    const double sn=arf_get_d(arb_midref(sine.value),ARF_RND_NEAR);
    if (max_abs(yc)>.71 || sc.lo<=0.0)
      throw std::runtime_error("fast threshold: point outside certified recurrence range");
    const double dy=up(std::max(std::abs(y-yc.lo),std::abs(yc.hi-y)));
    const double ds=up(std::max(std::abs(sn-sc.lo),std::abs(sc.hi-sn)));
    double previous=1.0,current=y,uprevious=0.0,ucurrent=1.0;
    double p=0.0,a=0.0,b=0.0;
    for (std::size_t k=0;k<g.size();++k) {
      const double j=static_cast<double>(k+1);
      if (std::abs(current)>2.0 || std::abs(ucurrent)>2.0)
        throw std::runtime_error("fast threshold: recurrence bootstrap failed");
      p += j*g[k]*current;
      a += j*j*g[k]*sn*ucurrent;
      b += j*j*j*g[k]*current;
      const double next=2.0*y*current-previous;
      const double unext=2.0*y*ucurrent-uprevious;
      previous=current;current=next;uprevious=ucurrent;ucurrent=unext;
    }
    if (!std::isfinite(p)||!std::isfinite(a)||!std::isfinite(b))
      throw std::runtime_error("fast threshold: nonfinite recurrence sum");
    const Interval u(std::ldexp(1.0,-53)),round_factor=Interval(4.0)*Interval(gamma);
    const Interval tiny=Interval(8.0)*Interval(static_cast<double>(g.size()))
      *Interval(std::numeric_limits<double>::denorm_min());
    // ep, ea, eb bound errors in P, -P', and -P''. The 1.43 bound for
    // 1/sin(theta) and the |cos(theta)|<=.71 guard cover the central domain.
    // tiny also covers subnormal rounding contributions.
    const double ep=(Interval(8.0)*u*Interval(s[1])+round_factor*Interval(s[0])
      +Interval(1.43)*Interval(s[1])*Interval(dy)+tiny*Interval(s[0])+tiny).hi;
    const double ea=(Interval(8.0)*u*Interval(s[2])+round_factor*Interval(s[1])
      +Interval(1.43)*Interval(s[1])*Interval(ds)
      +(Interval(3.0)*Interval(s[2])+Interval(3.0)*Interval(s[1]))*Interval(dy)
      +tiny*Interval(s[1])+tiny).hi;
    const double eb=(Interval(8.0)*u*Interval(s[3])+round_factor*Interval(s[2])
      +Interval(1.43)*Interval(s[3])*Interval(dy)+tiny*Interval(s[2])+tiny).hi;
    const Interval pp=centered(p,ep),aa=centered(a,ea),bb=centered(b,eb);
    const Interval y2=yc*yc;
    const Interval denominator=Interval(1.0)+Interval(4.0)*y2-Interval(8.0)*y2*y2;
    if (denominator.lo<.98)
      throw std::runtime_error("fast threshold: denominator lower bound failed");
    const Interval d1=Interval(8.0)*sc*yc*(Interval(4.0)*y2-Interval(1.0));
    const Interval d2=Interval(8.0)*(Interval(2.0)*y2-Interval(1.0))
      +Interval(16.0)*(Interval(8.0)*y2*y2-Interval(8.0)*y2+Interval(1.0));
    const Interval dd=denominator*denominator,ddd=dd*denominator;
    return {theta,-pp/denominator,
      aa/denominator+pp*d1/dd,
      bb/denominator-Interval(2.0)*aa*d1/dd+pp*d2/dd
        -Interval(2.0)*pp*d1*d1/ddd};
  }
};
// A theta subinterval and an upper bound for F everywhere on it.
struct Node {
  double left,right,upper;
  Point midpoint;
  bool operator<(const Node& other) const { return upper<other.upper; }
};
// Radius from a point to the farther endpoint, rounded outward.
inline double radius(double left,double right,double center) {
  return up(std::max(center-left,right-center));
}
// Taylor's theorem gives a global upper bound from midpoint derivatives.
// When F'' is negative throughout this subinterval, the concave quadratic
// bound F(p)+F'(p)^2/(-2*sup F'') is often substantially tighter.
inline double bound(Evaluator& evaluator,Point p,double left,double right) {
  const Interval r(radius(left,right,p.theta));
  const double ordinary=(Interval(p.value.hi)+Interval(max_abs(p.first))*r
    +Interval(.5)*Interval(evaluator.curvature)*r*r).hi;
  const double second_upper=(Interval(p.second.hi)+Interval(evaluator.third_bound)*r).hi;
  if (second_upper<0.0) {
    const Interval fp(max_abs(p.first));
    return std::min(ordinary,(Interval(p.value.hi)+fp*fp
      /(Interval(-2.0)*Interval(second_upper))).hi);
  }
  return ordinary;
}
// A bracket containing exactly one local maximum, established by F' signs
// and strictly negative curvature throughout the bracket.
struct RootCertificate {
  bool valid=false;
  double theta=0.0,left=0.0,right=0.0;
  Point point;
};
// Newton steps propose a location. Acceptance requires interval bounds,
// positive slope at the left endpoint, and negative slope at the right.
// Return valid=false when these conditions cannot be established.
inline RootCertificate stationary_certificate(Evaluator& evaluator,Point start,
    double domain_left,double domain_right,bool precise=false) {
  auto evaluate=[&](double theta) {
    return precise ? evaluator.precise(theta) : evaluator(theta);
  };
  Point p=precise ? evaluate(start.theta) : start;
  for (int step=0;step<12;++step) {
    const double curvature_mid=p.second.lo/2.0+p.second.hi/2.0;
    if (!(curvature_mid<0.0)) return {};
    const double first_mid=p.first.lo/2.0+p.first.hi/2.0;
    const double proposed=p.theta-first_mid/curvature_mid;
    if (!(proposed>domain_left && proposed<domain_right)) return {};
    if (proposed==p.theta) break;
    p=evaluate(proposed);
  }
  if (!(p.second.hi<0.0)) return {};
  double r=std::max(8.0*std::numeric_limits<double>::epsilon(),
    up(4.0*max_abs(p.first)/(-p.second.hi)));
  for (int attempt=0;attempt<12;++attempt,r*=2.0) {
    const double left=down(p.theta-r),right=up(p.theta+r);
    if (!(left>domain_left && right<domain_right)) return {};
    const double second_upper=(Interval(p.second.hi)
      +Interval(evaluator.third_bound)*Interval(radius(left,right,p.theta))).hi;
    if (!(second_upper<0.0)) return {};
    const Point a=evaluate(left),b=evaluate(right);
    if (a.first.lo>0.0 && b.first.hi<0.0)
      return {true,p.theta,left,right,p};
  }
  return {};
}
}  // namespace fast_detail

// Optional numerical diagnostics. Counts describe work, not sampler mixing.
// max_gradient_error_bound applies when the fast active-gradient path succeeds.
struct FastThresholdStats {
  long evaluations=0,nodes=0,precise_evaluations=0;
  bool stationary_certified=false,winner_separated=false,used_fallback=false;
  double stationary_left=0.0,stationary_right=0.0;
  double max_gradient_error_bound=0.0;
};

// Main entry point. g contains m original-mode slots with g[1]=0 (mode 2).
// Return the same maximum enclosure as the polynomial routine. A priority
// queue subdivides the interval with the largest remaining upper bound.
// With require_active_gradient=false only the value enclosure is populated.
// Otherwise, isolate and separate a winning location before returning its
// approximate envelope gradient. Ties or unresolved intervals can trigger
// the exact-polynomial fallback when allow_fallback=true.
inline ThresholdResult global_monotonicity_threshold_fast(
    const std::vector<double>& g,double left,double right,double half_width,
    double center=0.0,const ThresholdOptions& options=ThresholdOptions{},
    FastThresholdStats* stats=nullptr,bool allow_fallback=true,
    bool require_active_gradient=true) {
  using namespace fast_detail;
  try {
    if (g.size()<4 || g.size()>options.max_modes || options.max_modes>4096)
      throw std::invalid_argument("fast threshold: unsupported mode count");
    for(double x:g) if(!std::isfinite(x))
      throw std::invalid_argument("fast threshold: nonfinite coefficient");
    if(g[1]!=0.0) throw std::invalid_argument("fast threshold: second free coefficient must be zero");
    if(!std::isfinite(left)||!std::isfinite(right)||!std::isfinite(center)
      ||!std::isfinite(half_width)||!(left<right)||!(half_width>0.0)
      ||!std::isfinite(options.absolute_tolerance)||options.absolute_tolerance<=0.0
      ||options.relative_tolerance!=0.0 || options.max_precision_bits<96)
      throw std::invalid_argument("fast threshold: invalid domain or precision control");
    detail::Rational lr,rr,quarter,three_quarters;
    detail::endpoint_ratio(lr.value,left,center,half_width);
    detail::endpoint_ratio(rr.value,right,center,half_width);
    fmpq_set_si(quarter.value,1,4);fmpq_set_si(three_quarters.value,3,4);
    if(fmpq_cmp(lr.value,quarter.value)<0||fmpq_cmp(rr.value,three_quarters.value)>0)
      throw std::invalid_argument("fast threshold: padding must satisfy L>=2R");
    detail::Ball a,b,pi;
    arb_const_pi(pi.value,128);
    arb_set_fmpq(a.value,lr.value,128);arb_mul(a.value,a.value,pi.value,128);
    arb_set_fmpq(b.value,rr.value,128);arb_mul(b.value,b.value,pi.value,128);
    const Interval domain_a=ball_interval(a.value,128),domain_b=ball_interval(b.value,128);
    if(domain_a.hi>=domain_b.lo)
      throw std::runtime_error("fast threshold: no representable feasible theta point");
    Evaluator evaluator(g,options.max_precision_bits);
    // Feasible point values use inward-rounded endpoints for lower bounds.
    // The search covers outward-rounded endpoints so it omits no location.
    Point best=evaluator(domain_a.hi),right_point=evaluator(domain_b.lo);
    if(right_point.value.lo>best.value.lo) best=right_point;
    double lower=best.value.lo;
    std::priority_queue<Node> queue;
    auto insert=[&](double l,double r) {
      const double mid=l/2.0+r/2.0;
      if(!(mid>l&&mid<r)) throw std::runtime_error("fast threshold: subdivision exhausted");
      Point p=evaluator(mid);
      if(mid>=domain_a.hi&&mid<=domain_b.lo&&p.value.lo>lower) {lower=p.value.lo;best=p;}
      const double upper=bound(evaluator,p,l,r);
      queue.push({l,r,upper,p});
    };
    insert(domain_a.lo,domain_b.hi);
    long nodes=0;
    // This preliminary phase proposes a winner. The active-gradient path
    // subsequently separates the entire complement and checks the requested
    // final width. A value-only certificate needs the requested width here.
    const double preliminary_tolerance=require_active_gradient ?
      std::max(options.absolute_tolerance,1e-10) : options.absolute_tolerance;
    while(!queue.empty() &&
      (Interval(queue.top().upper)-Interval(lower)).hi>preliminary_tolerance) {
      const Node n=queue.top();queue.pop();
      if(n.upper<lower) continue;
      if(++nodes>200000) throw std::runtime_error("fast threshold: subdivision budget exhausted");
      const double mid=n.midpoint.theta;
      insert(n.left,mid);insert(mid,n.right);
      if(nodes%32==0 && best.second.hi<0.0) {
        const auto root=stationary_certificate(evaluator,best,domain_a.hi,domain_b.lo);
        if(root.valid && root.point.value.lo>lower) {best=root.point;lower=best.value.lo;}
      }
    }
    if(queue.empty()) throw std::runtime_error("fast threshold: lost global enclosure");
    if(!require_active_gradient) {
      ThresholdResult result;
      result.lower=lower;result.upper=std::max(lower,queue.top().upper);
      result.estimate=lower/2.0+result.upper/2.0;
      result.precision_bits=96;result.real_roots=-1;
      if(stats) {stats->evaluations=evaluator.evaluations;stats->nodes=nodes;
        stats->precise_evaluations=evaluator.precise_evaluations;}
      return result;
    }
    // Reserve most of the requested width for bounding the maximum between
    // points. Tighten point evaluations to one eighth of that tolerance.
    auto precise_value=[&](Point point) {
      while((Interval(point.value.hi)-Interval(point.value.lo)).hi>
        options.absolute_tolerance/8.0) {
        point=evaluator.precise(point.theta);
        if((Interval(point.value.hi)-Interval(point.value.lo)).hi<=
          options.absolute_tolerance/8.0) break;
        if(evaluator.precise_bits==options.max_precision_bits) break;
        evaluator.precise_bits=std::min(options.max_precision_bits,evaluator.precise_bits*2);
      }
      return point;
    };
    best=precise_value(best);
    lower=std::max(lower,best.value.lo);
    auto root=stationary_certificate(evaluator,best,domain_a.hi,domain_b.lo);
    const bool endpoint=best.theta==domain_a.hi || best.theta==domain_b.lo;
    const double jmax=static_cast<double>(g.size());
    // These absolute gradient tolerances govern location refinement. They
    // are numerical choices distinct from the threshold-value tolerance.
    const double gradient_lipschitz=(Interval(1.43)*Interval(jmax*jmax)/Interval(.98)
      +Interval(6.0)*Interval(jmax)/(Interval(.98)*Interval(.98))).hi;
    if(!endpoint && (!root.valid
        || up(gradient_lipschitz*radius(root.left,root.right,root.theta))>5e-7
        || (Interval(root.point.value.hi)-Interval(root.point.value.lo)).hi>
          options.absolute_tolerance/8.0)) {
      for(;;) {
        root=stationary_certificate(evaluator,best,domain_a.hi,domain_b.lo,true);
        if(root.valid
          && up(gradient_lipschitz*radius(root.left,root.right,root.theta))<=5e-7
          && (Interval(root.point.value.hi)-Interval(root.point.value.lo)).hi<=
            options.absolute_tolerance/8.0) break;
        if(evaluator.precise_bits==options.max_precision_bits) break;
        evaluator.precise_bits=std::min(options.max_precision_bits,evaluator.precise_bits*2);
      }
    }
    // A stationary certificate ensures the envelope gradient is evaluated at
    // a certified isolated critical point. Endpoints need no stationary test.
    if(!root.valid && !endpoint)
      throw std::runtime_error("fast threshold: winning stationary point unresolved");
    // Build a neighborhood with provable concavity for an interior winner,
    // or a one-sided monotone neighborhood for an endpoint winner.
    double winning_left,winning_right,upper;
    if(root.valid) {
      best=root.point;
      lower=std::max(lower,best.value.lo);
      const double r=(Interval(-best.second.hi)
        /(Interval(2.0)*Interval(evaluator.third_bound))).lo;
      winning_left=std::max(domain_a.lo,up(best.theta-r));
      winning_right=std::min(domain_b.hi,down(best.theta+r));
      if(!(winning_left<root.left && winning_right>root.right))
        throw std::runtime_error("fast threshold: winning curvature region unresolved");
      upper=bound(evaluator,best,winning_left,winning_right);
    } else {
      const bool at_left=best.theta==domain_a.hi;
      const double slope=at_left ? -best.first.hi : best.first.lo;
      if(!(slope>0.0))
        throw std::runtime_error("fast threshold: endpoint slope unresolved");
      const double r=(Interval(slope)/(Interval(2.0)*Interval(evaluator.curvature))).lo;
      winning_left=at_left ? domain_a.lo : std::max(domain_a.lo,up(best.theta-r));
      winning_right=at_left ? std::min(domain_b.hi,down(best.theta+r)) : domain_b.hi;
      if((Interval(slope)-Interval(evaluator.curvature)
          *Interval(radius(winning_left,winning_right,best.theta))).lo<=0.0)
        throw std::runtime_error("fast threshold: endpoint monotonicity unresolved");
      upper=precise_value(evaluator(at_left ? domain_a.lo : domain_b.hi)).value.hi;
    }
    const Point certified_winner=best;
    // Every interval outside the certified concave/monotone winning region
    // must have upper value strictly below a feasible point's lower value.
    // This is stronger than terminating at a small maximum-value gap.
    while(!queue.empty()) {
      const Node n=queue.top();queue.pop();
      if(n.upper<lower) continue;
      if(n.left>=winning_left && n.right<=winning_right) continue;
      if(++nodes>200000)
        throw std::runtime_error("fast threshold: competing maximum unresolved");
      if(n.left<winning_left && n.right>winning_left) {
        insert(n.left,winning_left);insert(winning_left,n.right);
      } else if(n.left<winning_right && n.right>winning_right) {
        insert(n.left,winning_right);insert(winning_right,n.right);
      } else {
        insert(n.left,n.midpoint.theta);insert(n.midpoint.theta,n.right);
      }
      if(lower>upper)
        throw std::runtime_error("fast threshold: competing point exceeds chosen winner");
    }
    if(!(upper>=lower) ||
      (Interval(upper)-Interval(lower)).hi>options.absolute_tolerance)
      throw std::runtime_error("fast threshold: final enclosure misses tolerance");
    best=certified_winner;
    detail::Ball input,cosine;
    arb_set_d(input.value,best.theta);arb_cos(cosine.value,input.value,96);
    ThresholdResult result;
    result.lower=lower;result.upper=upper;result.estimate=lower/2.0+upper/2.0;
    result.active_y=arf_get_d(arb_midref(cosine.value),ARF_RND_NEAR);
    result.gradient.resize(g.size());detail::envelope_gradient(result.gradient,result.active_y);
    double point_radius=0.0;
    if(root.valid) point_radius=radius(root.left,root.right,root.theta);
    else point_radius=up(std::max(domain_a.hi-domain_a.lo,domain_b.hi-domain_b.lo));
    point_radius=(Interval(point_radius)+Interval(2.0)*Interval(std::ldexp(1.0,-53))).hi;
    const double j=static_cast<double>(g.size());
    const double gradient_error=(
      (Interval(1.43)*Interval(j*j)/Interval(.98)
        +Interval(6.0)*Interval(j)/(Interval(.98)*Interval(.98)))
        *Interval(point_radius)
      +Interval(8.0)*Interval(std::ldexp(1.0,-53))*Interval(j*j)/Interval(.98)
      +Interval(128.0)*Interval(std::ldexp(1.0,-53))*Interval(j)
        /(Interval(.98)*Interval(.98))).hi;
    if(gradient_error>1e-6)
      throw std::runtime_error("fast threshold: envelope-gradient enclosure too wide");
    result.unique_active=true;
    result.precision_bits=evaluator.precise_evaluations ? evaluator.precise_bits : 96;
    result.real_roots=-1;
    if(stats) {stats->evaluations=evaluator.evaluations;stats->nodes=nodes;
      stats->precise_evaluations=evaluator.precise_evaluations;
      stats->stationary_certified=root.valid;stats->stationary_left=root.left;
      stats->stationary_right=root.right;stats->winner_separated=true;
      stats->max_gradient_error_bound=gradient_error;}
    return result;
  } catch(...) {
    // Polynomial root isolation handles flat objectives, active-location
    // ties, and precision limits of the fast path. Its value remains bounded,
    // but its gradient does not inherit the fast path's location-error bound.
    if(!allow_fallback) throw;
    if(stats) stats->used_fallback=true;
    auto result=global_monotonicity_threshold(g,left,right,half_width,center,options);
    if((Interval(result.upper)-Interval(result.lower)).hi>options.absolute_tolerance) {
      auto tighter=options;
      tighter.absolute_tolerance=options.absolute_tolerance/2.0;
      result=global_monotonicity_threshold(g,left,right,half_width,center,tighter);
      if((Interval(result.upper)-Interval(result.lower)).hi>options.absolute_tolerance)
        throw std::runtime_error("fast threshold: fallback misses strict width tolerance");
    }
    return result;
  }
}
}  // namespace hsgp_global
#endif
