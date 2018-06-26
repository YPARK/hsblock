#include <random>
#include "mathutil.hh"
#include "rcpp_util.hh"

#ifndef SAMPLER_HH_
#define SAMPLER_HH_

template <typename T>
struct discrete_sampler_t {
  using Index = typename T::Index;
  using Scalar = typename T::Scalar;

  explicit discrete_sampler_t(const Index k)
      : K(k), prob(k), unif01(zero_val, one_val), safe_exp(max_val) {
    prob.setZero();
  }

  template <typename Derived, typename RNG>
  Index operator()(const Eigen::MatrixBase<Derived>& log_prob, RNG& rng) {
    max_val = log_prob.maxCoeff();

    prob = log_prob.unaryExpr(safe_exp);
    prob /= prob.sum();

    Scalar u = unif01(rng);
    Scalar cum = 0.0;
    Index kk = 0;
    for (Index j = 0; j < prob.size(); ++j) {
      cum += prob(j);
      if (u <= cum) {
        kk = j;
        break;
      }
    }
    return kk;
  }

  const Index K;
  T prob;
  Scalar max_val;

  std::uniform_real_distribution<Scalar> unif01;

  struct safe_exp_op_t {
    explicit safe_exp_op_t(Scalar& _maxval) : maxval(_maxval) {}

    const Scalar operator()(const Scalar& x) const {
      return fasterexp(x - maxval);
    }

    Scalar& maxval;
  } safe_exp;

  constexpr static Scalar zero_val = 0.0;
  constexpr static Scalar one_val = 1.0;
};

#endif
