#ifndef HSB_DATA_HH_
#define HSB_DATA_HH_

#include "rcpp_util.hh"

struct tag_bernoulli {};
struct tag_poisson {};
struct tag_gaussian {};

template <typename T>
struct data_traits {
  using Distrib = typename T::Distrib;  // distribution tag
  using Unit = typename T::Unit;        // unit data type
};

template <typename T>
using distrib_tag = typename data_traits<T>::Distrib;

#include "hsb_bernoulli.hh"

template <typename DT>
void clear(DT& dt) {
  impl_clear(dt, distrib_tag<DT>());
}

template <typename DT>
void dump(DT& dt) {
  impl_dump(dt, distrib_tag<DT>());
}

template <typename DT>
void eval_delta_score(DT& dt) {
  impl_eval_delta_score(dt, distrib_tag<DT>());
}

template <typename DT>
void merge_left_right_delta(const DT& left, const DT& right, DT& intern) {
  impl_merge_left_right_delta(left, right, intern, distrib_tag<DT>());
}

template <typename DT>
void eval_score(DT& dt) {
  impl_eval_score(dt, distrib_tag<DT>());
}

template <typename DT>
void resolve_delta(DT& dt) {
  impl_resolve_delta(dt, distrib_tag<DT>());
}

template <typename DT, typename Scalar>
void update_delta(DT& dt, const Scalar rate) {
  impl_update_delta(dt, rate, distrib_tag<DT>());
}

#endif
