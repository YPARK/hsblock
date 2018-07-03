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
#include "hsb_poisson.hh"

template <typename DT>
void clear(DT& dt) {
  impl_clear(dt, distrib_tag<DT>());
}

template <typename DT, typename Scalar>
void init_param_intern(DT& dt, DT& dt_left, DT& dt_right, const Scalar level) {
  impl_init_param_intern(dt, dt_left, dt_right, level, distrib_tag<DT>());
}

template <typename DT, typename Scalar>
void init_param(DT& dt, const Scalar level) {
  impl_init_param(dt, level, distrib_tag<DT>());
}

template <typename DT>
void init_param_null(DT& dt) {
  impl_init_param_null(dt, distrib_tag<DT>());
}

template <typename DT>
void dump(DT& dt) {
  impl_dump(dt, distrib_tag<DT>());
}

////////////////////////////
// variational parameters //
////////////////////////////

template <typename DT, typename Scalar>
void update(DT& dt, const Scalar rate) {
  impl_update(dt, rate, distrib_tag<DT>());
}

////////////////////////////////
// evaluation of local scores //
////////////////////////////////

template <typename DT>
void eval_delta_score(DT& dt) {
  impl_eval_delta_score(dt, distrib_tag<DT>());
}

template <typename DT>
void eval_delta_left(DT& dt, const DT& right) {
  impl_eval_delta_left(dt, right, distrib_tag<DT>());
}

template <typename DT>
void eval_delta_right(DT& dt, const DT& left) {
  impl_eval_delta_right(dt, left, distrib_tag<DT>());
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
void eval_stat_score(DT& dt) {
  impl_eval_stat_score(dt, distrib_tag<DT>());
}

#endif
