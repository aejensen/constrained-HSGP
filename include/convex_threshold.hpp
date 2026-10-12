#ifndef CONVEX_HSGP_THRESHOLD_HPP
#define CONVEX_HSGP_THRESHOLD_HPP

// Convexity threshold for f=a0+sum_j a_j*phi_j with gamma=-a1.
// The input vector has m entries: g[0]=a0, g[j-1]=a_j for j=2,...,m.
// For w=cos(pi*(x+L)/(2L)), the event f''(x)>=0 is equivalent to
//   gamma >= max_w Q(w),  Q(w)=sum_{j=2}^m j^2*g[j-1]*U_{j-1}(w).
// R centers the interval before calling this routine, so c=0 here. Require
// 2<=m<=32 and -L<left<right<L. The intercept never enters the threshold.
//
// The usual path bounds Q using Bernstein coefficients and subdivision.
// A rational-polynomial/Arb fallback isolates all stationary candidates.
// upper-lower is at most the requested tolerance, in coefficient units.
// gradient is an approximate active-location envelope vector, not the
// derivative of the rounded upper bound. Compile without fast-math or FMA
// contraction so outward rounding is preserved.

#include <flint/arb.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/arb_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <gmp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

namespace convex_hsgp {
// Outward rounding for interval additions and products in double precision.
inline double down(double x) { return std::nextafter(x, -INFINITY); }
inline double up(double x) { return std::nextafter(x, INFINITY); }
struct Interval { double lo, hi; Interval(double x = 0):lo(x),hi(x) {}
  Interval(double a,double b):lo(a),hi(b) {} };
inline Interval operator+(Interval a,Interval b) {return {down(a.lo+b.lo),up(a.hi+b.hi)};}
inline Interval operator*(Interval a,Interval b) {
  double p[]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
  return {down(*std::min_element(p,p+4)),up(*std::max_element(p,p+4))};
}
// Own an Arb ball, including its initialization, copying, and cleanup.
struct Ball {
  arb_t v;
  Ball() {arb_init(v);}
  Ball(const Ball& b) {arb_init(v);arb_set(v,b.v);}
  Ball& operator=(const Ball& b) {arb_set(v,b.v);return *this;}
  ~Ball() {arb_clear(v);}
};
// Convert an Arb ball to outward-rounded double endpoints.
inline Interval enclosure(const Ball& b) {
  arf_t lo,hi; arf_init(lo);arf_init(hi);
  arb_get_interval_arf(lo,hi,b.v,128);
  Interval ans(arf_get_d(lo,ARF_RND_FLOOR),arf_get_d(hi,ARF_RND_CEIL));
  arf_clear(lo);arf_clear(hi); return ans;
}
// Binomial coefficients used in the power-to-Bernstein transformation.
inline std::uint64_t choose(int n,int k) {
  std::uint64_t ans=1;
  for(int i=1;i<=k;++i) ans=ans*(n-i+1)/i;
  return ans;
}
// active is in cosine coordinates. gradient has the input g ordering, with
// a zero intercept entry. nodes=-1 indicates use of the root fallback.
struct Result {double lower,upper,active;std::vector<double> gradient;int nodes;};
inline Result root_certificate(const std::vector<double>& g,double left,double right,double L,double tolerance);
// Precompute, for fixed geometry, a map from g to Bernstein coefficients.
// Each coefficient map is enclosed with Arb before conversion to doubles.
struct Oracle {
  int m,n; double left,right,L,A,B;
  std::vector<std::vector<Interval>> maps;
  Oracle(int modes,double xl,double xr,double pad):m(modes),n(modes-1),left(xl),right(xr),L(pad) {
    if(m<2 || m>32 || !(left<right) || !(left>-L && right<L))
      throw std::domain_error("Convex oracle requires 2 <= rank <= 32 and an interior interval.");
    Ball pi,bl,theta,a,b,width;
    arb_const_pi(pi.v,128);arb_set_d(bl.v,L);
    auto coordinate=[&](double x,Ball& out) {
      arb_set_d(theta.v,x);arb_add(theta.v,theta.v,bl.v,128);
      arb_mul(theta.v,theta.v,pi.v,128);arb_div(theta.v,theta.v,bl.v,128);
      arb_mul_2exp_si(theta.v,theta.v,-1);arb_cos(out.v,theta.v,128);
    };
    coordinate(right,a);coordinate(left,b);arb_sub(width.v,b.v,a.v,128);
    A=enclosure(a).lo;B=enclosure(b).hi;
    std::vector<Ball> ap(m),wp(m);
    arb_one(ap[0].v);arb_one(wp[0].v);
    for(int i=1;i<m;++i) {arb_mul(ap[i].v,ap[i-1].v,a.v,128);arb_mul(wp[i].v,wp[i-1].v,width.v,128);}
    // Generate U_j by U_j(w)=2*w*U_{j-1}(w)-U_{j-2}(w), starting at 1,2w.
    // Integer coefficients fit in the chosen integer types for m<=32.
    std::vector<std::vector<std::int64_t>> U(m,std::vector<std::int64_t>(m,0));
    U[0][0]=1;U[1][1]=2;
    for(int j=2;j<m;++j) for(int i=0;i<=j;++i)
      U[j][i]=(i?2*U[j-1][i-1]:0)-U[j-2][i];
    // Substitute w=a+(b-a)*t, t in [0,1], then convert powers of t to
    // degree-n Bernstein coefficients. Their convex-hull property bounds Q
    // over the full interval by the minimum/maximum coefficient bounds.
    maps.resize(m-1,std::vector<Interval>(m));
    for(int j=1;j<m;++j) {
      std::vector<Ball> power(m);Ball term,sum;
      for(int k=0;k<=j;++k) for(int i=k;i<=j;++i) if(U[j][i]) {
        arb_mul(term.v,ap[i-k].v,wp[k].v,128);
        arb_mul_si(term.v,term.v,U[j][i],128);
        arb_mul_ui(term.v,term.v,choose(i,k),128);
        arb_add(power[k].v,power[k].v,term.v,128);
      }
      for(int l=0;l<m;++l) {
        arb_zero(sum.v);
        for(int k=0;k<=std::min(l,j);++k) {
          arb_mul_ui(term.v,power[k].v,choose(l,k),128);
          arb_div_ui(term.v,term.v,choose(n,k),128);
          arb_add(sum.v,sum.v,term.v,128);
        }
        arb_mul_ui(sum.v,sum.v,(j+1)*(j+1),128);
        maps[j-1][l]=enclosure(sum);
      }
    }
  }
// Ordinary recurrence evaluation used only to refine the active location.
  double value(const std::vector<double>& g,double w) const {
    double u0=1,u1=2*w,ans=0;
    for(int j=2;j<=m;++j) {
      double u=(j==2?u1:2*w*u1-u0);
      if(j>2) {u0=u1;u1=u;}
      ans+=j*j*g[j-1]*u;
    }
    return ans;
  }
  // Bound the maximum, splitting the interval with the highest upper bound.
  // Subdivision uses the de Casteljau averaging scheme at t=1/2.
  Result operator()(const std::vector<double>& g,double tolerance) const {
    if(static_cast<int>(g.size())!=m || !(tolerance>0)) throw std::domain_error("Oracle input dimension/tolerance.");
    struct Node {std::vector<Interval> cf;double lo,hi,upper;};
    auto make_node=[](std::vector<Interval> cf,double lo,double hi) {
      double ub=-INFINITY;for(auto x:cf) ub=std::max(ub,x.hi);
      return Node{std::move(cf),lo,hi,ub};
    };
    auto compare=[](const Node& a,const Node& b) {return a.upper<b.upper;};
    std::priority_queue<Node,std::vector<Node>,decltype(compare)> queue(compare);
    std::vector<Interval> cf(m);
    for(int j=2;j<=m;++j) {
      if(!std::isfinite(g[j-1])) throw std::domain_error("Nonfinite residual coefficient.");
      for(int k=0;k<m;++k) cf[k]=cf[k]+Interval(g[j-1])*maps[j-2][k];
    }
    // Endpoint Bernstein coefficients equal endpoint function values, so
    // their lower endpoints supply feasible lower bounds for the maximum.
    double lower=cf.front().lo,active=A;
    if(cf.back().lo>lower) {lower=cf.back().lo;active=B;}
    queue.push(make_node(cf,A,B));int nodes=0;
    while(queue.top().upper-lower>tolerance) {
      // Switch algorithms if repeated subdivision is no longer economical.
      if(++nodes>256) return root_certificate(g,left,right,L,tolerance);
      Node top=queue.top();queue.pop();
      std::vector<Interval> work=top.cf,lc(m),rc(m);
      lc[0]=work[0];rc[n]=work[n];
      for(int d=1;d<m;++d) {
        for(int k=0;k<m-d;++k) work[k]=(work[k]+work[k+1])*Interval(.5);
        lc[d]=work[0];rc[n-d]=work[n-d];
      }
      double mid=(top.lo+top.hi)/2;
      if(lc[n].lo>lower) {lower=lc[n].lo;active=mid;}
      queue.push(make_node(std::move(lc),top.lo,mid));
      queue.push(make_node(std::move(rc),mid,top.hi));
    }
    double upper=queue.top().upper;
    if(!std::isfinite(lower) || !std::isfinite(upper)) throw std::domain_error("Nonfinite maximum enclosure.");
    // Refine the active point for the envelope derivative. The enclosure
    // remains independent of this ordinary floating-point refinement.
    if(active>A && active<B) {
      double w=active;
      for(int iter=0;iter<12;++iter) {
        double u0=1,u1=2*w,d0=0,d1=2,e0=0,e1=0,first=0,second=0;
        for(int j=2;j<=m;++j) {
          double u=u1,d=d1,e=e1;
          if(j>2) {u=2*w*u1-u0;d=2*u1+2*w*d1-d0;e=4*d1+2*w*e1-e0;
            u0=u1;u1=u;d0=d1;d1=d;e0=e1;e1=e;}
          first+=j*j*g[j-1]*d;second+=j*j*g[j-1]*e;
        }
        if(!(second<0)) break;
        double next=w-first/second;
        if(!(next>A && next<B) || value(g,next)<value(g,w)-tolerance) break;
        w=next;
      }
      if(value(g,w)>=value(g,active)-tolerance) active=w;
    }
    // dQ/dg_j=j^2*U_{j-1}(active). At tied maxima, this chooses one active
    // vector rather than asserting differentiability of the envelope.
    std::vector<double> gradient(m,0);double u0=1,u1=2*active;
    for(int j=2;j<=m;++j) {
      double u=j==2?u1:2*active*u1-u0;
      if(j>2) {u0=u1;u1=u;}
      gradient[j-1]=j*j*u;
    }
    return {lower,upper,active,gradient,nodes};
  }
};
// Build Q exactly from the binary-rational inputs, remove repeated factors
// from Q', and enclose all its real roots with Arb. Include both endpoints
// and any root whose membership is uncertain. Such an uncertain root may
// increase the upper bound, but cannot increase the feasible lower bound.
// Constant polynomials are covered by endpoint evaluation alone.
inline Result root_certificate(const std::vector<double>& g,double left,double right,double L,double tolerance) {
  fmpq_poly_t q,u0,u1,next,term,dq;
  fmpq_poly_init(q);fmpq_poly_init(u0);fmpq_poly_init(u1);
  fmpq_poly_init(next);fmpq_poly_init(term);fmpq_poly_init(dq);
  fmpq_t scalar,xl,xr,pad,ratio;
  fmpq_init(scalar);fmpq_init(xl);fmpq_init(xr);fmpq_init(pad);fmpq_init(ratio);
  mpq_t tmp;mpq_init(tmp);
  auto set_double=[&](fmpq_t dest,double x) {mpq_set_d(tmp,x);fmpq_set_mpq(dest,tmp);};
  set_double(xl,left);set_double(xr,right);set_double(pad,L);
  fmpq_poly_set_coeff_si(u0,0,1);fmpq_poly_set_coeff_si(u1,1,2);
  for(int j=2;j<=static_cast<int>(g.size());++j) {
    set_double(scalar,g[j-1]);fmpq_mul_ui(scalar,scalar,j*j);
    fmpq_poly_scalar_mul_fmpq(term,u1,scalar);fmpq_poly_add(q,q,term);
    fmpq_poly_shift_left(next,u1,1);fmpq_poly_scalar_mul_si(next,next,2);
    fmpq_poly_sub(next,next,u0);fmpq_poly_set(u0,u1);fmpq_poly_set(u1,next);
  }
  fmpq_poly_derivative(dq,q);
  fmpz_poly_t p,dp,gcd,sf;
  fmpz_poly_init(p);fmpz_poly_init(dp);fmpz_poly_init(gcd);fmpz_poly_init(sf);
  fmpq_poly_get_numerator(p,dq);
  if(fmpz_poly_degree(p)>0) {
    fmpz_poly_derivative(dp,p);fmpz_poly_gcd(gcd,p,dp);
    fmpz_poly_divexact(sf,p,gcd);fmpz_poly_primitive_part(p,sf);
  }
  int degree=std::max<slong>(0,fmpz_poly_degree(p));
  Result result;bool resolved=false;
  // Raise precision until the maximum-value enclosure meets tolerance.
  // Failure at the supported precision limit is reported to the caller.
  for(long bits=128;bits<=1024;bits*=2) {
    Ball a,b;fmpq_add(ratio,xr,pad);fmpq_div(ratio,ratio,pad);fmpq_div_2exp(ratio,ratio,1);
    arb_cos_pi_fmpq(a.v,ratio,bits);
    fmpq_add(ratio,xl,pad);fmpq_div(ratio,ratio,pad);fmpq_div_2exp(ratio,ratio,1);
    arb_cos_pi_fmpq(b.v,ratio,bits);
    arb_ptr roots=_arb_vec_init(degree);
    slong nr=degree?arb_fmpz_poly_real_roots(roots,p,0,bits):0;
    std::vector<Ball> points;std::vector<bool> inside;
    points.push_back(a);points.push_back(b);inside.push_back(true);inside.push_back(true);
    auto aa=enclosure(a),bb=enclosure(b);
    for(slong i=0;i<nr;++i) {
      Ball root;arb_set(root.v,roots+i);auto rr=enclosure(root);
      if(rr.hi<aa.lo || rr.lo>bb.hi) continue;
      points.push_back(root);inside.push_back(rr.lo>aa.hi && rr.hi<bb.lo);
    }
    arb_poly_t aq;arb_poly_init(aq);arb_poly_set_fmpq_poly(aq,q,bits);
    double lo=-INFINITY,hi=-INFINITY,best=-INFINITY,active=0;
    for(std::size_t i=0;i<points.size();++i) {
      Ball val;arb_poly_evaluate(val.v,aq,points[i].v,bits);auto vv=enclosure(val);
      if(inside[i]) lo=std::max(lo,vv.lo);
      hi=std::max(hi,vv.hi);
      double mid=arf_get_d(arb_midref(val.v),ARF_RND_NEAR);
      if(mid>best) {best=mid;active=arf_get_d(arb_midref(points[i].v),ARF_RND_NEAR);}
    }
    arb_poly_clear(aq);_arb_vec_clear(roots,degree);
    if(std::isfinite(lo) && std::isfinite(hi) && hi-lo<=tolerance) {
      result={lo,hi,active,std::vector<double>(g.size(),0),-1};resolved=true;break;
    }
  }
  fmpq_poly_clear(q);fmpq_poly_clear(u0);fmpq_poly_clear(u1);fmpq_poly_clear(next);
  fmpq_poly_clear(term);fmpq_poly_clear(dq);fmpq_clear(scalar);fmpq_clear(xl);
  fmpq_clear(xr);fmpq_clear(pad);fmpq_clear(ratio);mpq_clear(tmp);
  fmpz_poly_clear(p);fmpz_poly_clear(dp);fmpz_poly_clear(gcd);fmpz_poly_clear(sf);
  if(!resolved) throw std::domain_error("Polynomial maximum unresolved at 1024 bits.");
  double u0d=1,u1d=2*result.active;
  for(int j=2;j<=static_cast<int>(g.size());++j) {
    double u=j==2?u1d:2*result.active*u1d-u0d;
    if(j>2) {u0d=u1d;u1d=u;}
    result.gradient[j-1]=j*j*u;
  }
  return result;
}
// Reuse geometry maps across Stan evaluations on the same worker thread.
// The cache depends only on rank and interval geometry, never on g or data.
inline const Oracle& oracle(int m,double left,double right,double L) {
  static thread_local std::vector<Oracle> cache;
  for(const auto& o:cache) if(o.m==m && o.left==left && o.right==right && o.L==L) return o;
  cache.emplace_back(m,left,right,L);return cache.back();
}
}
#endif
