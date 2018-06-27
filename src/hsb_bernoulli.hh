#ifndef HSB_BERNOULLI_HH_
#define HSB_BERNOULLI_HH_

#include "mathutil.hh"
#include "rcpp_util.hh"

template <typename T>
struct hsb_bern_t {
  using Unit = T;
  using Distrib = tag_bernoulli;

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
  T var_hole;  // variational parameter for holes
};

template <typename DT>
inline void impl_clear(DT& data, tag_bernoulli) {
  data.a0 = 1.0;
  data.b0 = 1.0;
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

template <typename DT>
void impl_dump(DT& data, tag_bernoulli) {
  Rcpp::Rcerr << "E: " << data.stat_edge;
  Rcpp::Rcerr << " T: " << data.stat_total;
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
                                        DT& intern, tag_bernoulli) {
  intern.delta_stat_dik = left.delta_stat_dik + right.delta_stat_dik;
  intern.delta_stat_nk = left.delta_stat_nk + right.delta_stat_nk;
  intern.delta_stat_nik = left.delta_stat_nik + right.delta_stat_nik;
}

///////////////////////////////////////////
// Gam(dE + A) Gam(dH + B) Gam(A + B)	 //
// ----------------------- ------------- //
// Gam(dE + dH + A + B)    Gam(A) Gam(B) //
///////////////////////////////////////////
template <typename DT>
inline void impl_eval_delta_score(DT& D, tag_bernoulli) {
  static typename DT::Unit TOL = 1e-4;

#ifdef DEBUG
  ASSERT(D.delta_stat_dik <= (D.delta_stat_nk + TOL),
         "dik <= nk : " << D.delta_stat_dik << " vs " << D.delta_stat_nk);
  ASSERT(D.stat_edge <= (D.stat_total + TOL),
         "[bottom] edge <= total : " << D.stat_edge << " vs " << D.stat_total);
#endif

  D.var_edge = D.stat_edge + D.a0 + TOL;
  D.var_hole = D.stat_total - D.stat_edge + D.b0 + TOL;

  // First term
  D.delta_score = fasterlgamma(D.delta_stat_dik + D.var_edge);
  D.delta_score +=
      fasterlgamma(D.delta_stat_nk - D.delta_stat_dik + D.var_hole);
  D.delta_score -= fasterlgamma(D.delta_stat_nk + D.var_edge + D.var_hole);

  // Second term
  D.delta_score += fasterlgamma(D.var_edge + D.var_hole);
  D.delta_score -= fasterlgamma(D.var_edge);
  D.delta_score -= fasterlgamma(D.var_hole);
}

template <typename DT>
inline void impl_eval_delta_left(DT& D, const DT& right, tag_bernoulli) {
  static typename DT::Unit TOL = 1e-4;

#ifdef DEBUG
  ASSERT(right.delta_stat_dik <= (TOL + right.delta_stat_nk),
         "right.dik <= nk");
  ASSERT(D.stat_edge <= (TOL + D.stat_total),
         "[left] edge <= total : " << D.stat_edge << " vs " << D.stat_total);
#endif

  D.var_edge = D.stat_edge + D.a0 + TOL;
  D.var_hole = D.stat_total - D.stat_edge + D.b0 + TOL;

  // First term
  D.delta_score_left = fasterlgamma(right.delta_stat_dik + D.var_edge);
  D.delta_score_left +=
      fasterlgamma(right.delta_stat_nk - right.delta_stat_dik + D.var_hole);
  D.delta_score_left -=
      fasterlgamma(right.delta_stat_nk + D.var_edge + D.var_hole);

  // Second term
  D.delta_score_left += fasterlgamma(D.var_edge + D.var_hole);
  D.delta_score_left -= fasterlgamma(D.var_edge);
  D.delta_score_left -= fasterlgamma(D.var_hole);
}

template <typename DT>
inline void impl_eval_delta_right(DT& D, const DT& left, tag_bernoulli) {
  static typename DT::Unit TOL = 1e-4;

#ifdef DEBUG
  ASSERT(left.delta_stat_dik <= (TOL + left.delta_stat_nk), "left.dik <= nk");
  ASSERT(D.stat_edge <= (TOL + D.stat_total),
         "[right] edge <= total : " << D.stat_edge << " vs " << D.stat_total);
#endif

  D.var_edge = D.stat_edge + D.a0 + TOL;
  D.var_hole = D.stat_total - D.stat_edge + D.b0 + TOL;

  // First term
  D.delta_score_right = fasterlgamma(left.delta_stat_dik + D.var_edge);
  D.delta_score_right +=
      fasterlgamma(left.delta_stat_nk - left.delta_stat_dik + D.var_hole);
  D.delta_score_right -=
      fasterlgamma(left.delta_stat_nk + D.var_edge + D.var_hole);

  // Second term
  D.delta_score_right += fasterlgamma(D.var_edge + D.var_hole);
  D.delta_score_right -= fasterlgamma(D.var_edge);
  D.delta_score_right -= fasterlgamma(D.var_hole);
}

///////////////////
// Gam(A + B)	 //
// ------------- //
// Gam(A) Gam(B) //
///////////////////
template <typename DT>
inline void impl_eval_score(DT& D, tag_bernoulli) {
  static typename DT::Unit TOL = 1e-4;
  D.var_edge = D.stat_edge + D.a0 + TOL;
  D.var_hole = D.stat_total - D.stat_edge + D.b0 + TOL;

  D.score = fasterlgamma(D.var_edge);
  D.score += fasterlgamma(D.var_hole);
  D.score -= fasterlgamma(D.var_edge + D.var_hole);
}

#endif
