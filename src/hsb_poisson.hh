#ifndef HSB_POISSON_HH_
#define HSB_POISSON_HH_

#include "mathutil.hh"
#include "rcpp_util.hh"

template <typename T>
struct hsb_pois_t {
  using Unit = T;
  using Distrib = tag_poisson;

  T stat_edge;   // stat for edge
  T stat_total;  // stat for total
  T score;       // log-likelihood tyoe of score

  T delta_stat_dik;  // delta stat for edge
  T delta_stat_nik;  // delta stat for total
  T delta_stat_nk;   // delta stat for total

  T delta_score;        // delta score
  T delta_score_left;   // delta score toward the left
  T delta_score_right;  // delta score toward the right

  T a0;        // prior for edges
  T b0;        // prior for holes
  T var_edge;  // variational parameter for edges
  T var_tot;   // variational parameter for holes
};

template <typename DT>
inline void impl_clear(DT& data, tag_poisson) {
  data.stat_edge = 0.0;
  data.stat_total = 0.0;
  data.score = 0.0;
  data.delta_stat_dik = 0.0;
  data.delta_stat_nik = 0.0;
  data.delta_stat_nk = 0.0;
  data.delta_score = 0.0;
  data.delta_score_left = 0.0;
  data.delta_score_right = 0.0;
}

template <typename DT, typename Scalar>
inline void impl_init_param(DT& data, const Scalar level, tag_poisson) {
  using Unit = typename DT::Unit;

  const Unit scale = 5.0;
  const Unit diminishing = 2.0;
  const Unit min_nk = 10.0;
  const Unit TOL = 1e-4;
  const Unit one_val = 1.0;
  const Unit half_val = 0.5;

  data.a0 = one_val;
  data.b0 = TOL;

  const Unit pr = std::pow(diminishing, -level);

  const Unit Nk = std::max(data.delta_stat_nk, min_nk);
  const Unit Tot = Nk * (Nk - one_val) * half_val;

  data.var_edge = one_val + Tot * pr * scale;
  data.var_tot = one_val + Tot;
}

template <typename DT, typename Scalar>
inline void impl_init_param_intern(DT& data, DT& data_left, DT& data_right,
                                   const Scalar level, tag_poisson) {
  using Unit = typename DT::Unit;

  const Unit scale = 5.0;
  const Unit diminishing = 2.0;
  const Unit TOL = 1e-4;
  const Unit one_val = 1.0;

  data.a0 = one_val;
  data.b0 = TOL;

  const Unit pr = std::pow(diminishing, -level);

  const Unit nL = data_left.delta_stat_nk;
  const Unit nR = data_right.delta_stat_nk;
  const Unit Tot = nL * nR;

  data.var_edge = one_val + Tot * pr * scale;
  data.var_tot = one_val + Tot;
}

template <typename DT>
inline void impl_init_param_null(DT& data, tag_poisson) {
  using Unit = typename DT::Unit;
  const Unit one_val = 1.0;

  data.a0 = one_val;
  data.b0 = one_val + one_val;
  data.var_edge = data.stat_edge;
  data.var_tot = data.stat_total;
}

template <typename DT>
void impl_dump(DT& data, tag_poisson) {
  Rcpp::Rcerr << "E: " << data.stat_edge;
  Rcpp::Rcerr << " T: " << data.stat_total;
  Rcpp::Rcerr << " var_edge: " << data.var_edge;
  Rcpp::Rcerr << " var_tot: " << data.var_tot;
  Rcpp::Rcerr << " d_ik: " << data.delta_stat_dik;
  Rcpp::Rcerr << " n_ik: " << data.delta_stat_nik;
  Rcpp::Rcerr << " n_k: " << data.delta_stat_nk;
  Rcpp::Rcerr << " dS: " << data.delta_score;
  Rcpp::Rcerr << " dSL: " << data.delta_score_left;
  Rcpp::Rcerr << " dSR: " << data.delta_score_right;
  Rcpp::Rcerr << std::endl;
}

template <typename DT>
inline void impl_merge_left_right_delta(const DT& left, const DT& right,
                                        DT& intern, tag_poisson) {
  intern.delta_stat_dik = left.delta_stat_dik + right.delta_stat_dik;
  intern.delta_stat_nik = left.delta_stat_nik + right.delta_stat_nik;
  intern.delta_stat_nk = left.delta_stat_nk + right.delta_stat_nk;
}

// Gamma(dE + A)      B^A
// ----------------- ---------
// (dT + B)^(dE + A)  Gamma(A)
template <typename DT>
inline void impl_eval_delta_score(DT& D, tag_poisson) {
  static typename DT::Unit TOL = 1e-4;

  // Gamma terms
  D.delta_score = fasterlgamma(D.delta_stat_dik + D.var_edge + TOL);
  D.delta_score -= fasterlgamma(D.var_edge + TOL);

  // power terms
  D.delta_score += fasterlog(D.var_tot + TOL) * (D.var_edge + TOL);
  D.delta_score -= fasterlog(D.var_tot + D.delta_stat_nk + TOL) *
                   (D.var_edge + D.delta_stat_dik + TOL);

  // dump(D);
}

template <typename DT>
inline void impl_eval_delta_left(DT& D, const DT& right, tag_poisson) {
  static typename DT::Unit TOL = 1e-4;

  // Gamma terms
  D.delta_score_left = fasterlgamma(right.delta_stat_dik + D.var_edge + TOL);
  D.delta_score_left -= fasterlgamma(D.var_edge + TOL);

  // power terms
  D.delta_score_left += fasterlog(D.var_tot + TOL) * (D.var_edge + TOL);
  D.delta_score_left -= fasterlog(D.var_tot + right.delta_stat_nk + TOL) *
                        (D.var_edge + right.delta_stat_dik + TOL);
}

template <typename DT>
inline void impl_eval_delta_right(DT& D, const DT& left, tag_poisson) {
  static typename DT::Unit TOL = 1e-4;

  // Gamma terms
  D.delta_score_right = fasterlgamma(left.delta_stat_dik + D.var_edge + TOL);
  D.delta_score_right -= fasterlgamma(D.var_edge + TOL);

  // power terms
  D.delta_score_right += fasterlog(D.var_tot + TOL) * (D.var_edge + TOL);
  D.delta_score_right -= fasterlog(D.var_tot + left.delta_stat_nk + TOL) *
                         (D.var_edge + left.delta_stat_dik + TOL);
}

template <typename DT, typename Scalar>
inline void impl_update(DT& D, const Scalar rate, tag_poisson) {
  D.var_edge = D.var_edge * (1.0 - rate);
  D.var_edge += rate * (D.stat_edge + D.a0);
  D.var_tot = D.var_tot * (1.0 - rate);
  D.var_tot += rate * (D.stat_total + D.b0);
}

////////////////////////////////////
// b0^(a0)    Gam(a0 + Edge)      //
// ------- ---------------------- //
// Gam(a0) (b0 + Tot)^(a0 + Edge) //
////////////////////////////////////

template <typename DT>
inline void impl_eval_score(DT& D, tag_poisson) {
  static typename DT::Unit TOL = 1e-4;

  D.score = D.a0 * fasterlog(D.b0 + TOL);
  D.score -= fasterlgamma(D.a0 + TOL);
  D.score += fasterlgamma(D.a0 + D.var_edge + TOL);
  D.score -= (D.a0 + D.var_edge) * fasterlog(D.b0 + D.var_tot + TOL);
}

template <typename DT>
inline void impl_eval_stat_score(DT& D, tag_poisson) {
  static typename DT::Unit TOL = 1e-4;

  D.score = D.a0 * fasterlog(D.b0 + TOL);
  D.score -= fasterlgamma(D.a0 + TOL);
  D.score += fasterlgamma(D.a0 + D.stat_edge + TOL);
  D.score -= (D.a0 + D.stat_edge) * fasterlog(D.b0 + D.stat_total + TOL);
}

#endif
