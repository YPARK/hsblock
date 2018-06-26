#ifndef HSB_BERNOULLI_HH_
#define HSB_BERNOULLI_HH_

#include "rcpp_util.hh"

template <typename T>
struct hsb_bern_t {
  using Unit = T;
  using Distrib = tag_bernoulli;

  T stat_edge;          // stat for edge
  T stat_total;         // stat for total
  T score;              // log-likelihood tyoe of score
  T delta_stat_edge;    // delta stat for edge
  T delta_stat_total;   // delta stat for total
  T delta_score;        // delta score
  T delta_score_left;   // delta score toward the left
  T delta_score_right;  // delta score toward the right

  T a0;        // prior for edges
  T b0;        // prior for holes
  T var_edge;  // variational parameter for edges
  T var_hole;  // variational parameter for holes
};

template <typename DT>
void impl_clear(DT& data, tag_bernoulli) {
  data.stat_edge = 0.0;
  data.stat_total = 0.0;
  data.score = 0.0;
  data.delta_stat_edge = 0.0;
  data.delta_stat_total = 0.0;
  data.delta_score_left = 0.0;
  data.delta_score_right = 0.0;
}

template <typename DT>
void impl_dump(DT& data, tag_bernoulli) {
  Rcpp::Rcerr << "E: " << data.stat_edge;
  Rcpp::Rcerr << " T: " << data.stat_total;
  Rcpp::Rcerr << " dE: " << data.delta_stat_edge;
  Rcpp::Rcerr << " dT: " << data.delta_stat_total;
  Rcpp::Rcerr << " score: " << data.score;
  Rcpp::Rcerr << " scoreL: " << data.delta_score_left;
  Rcpp::Rcerr << " scoreR: " << data.delta_score_right;
  Rcpp::Rcerr << std::endl;
}

// Delta statistics to pass over
template<typename T>
struct hsb_bern_delta_t {
  explicit hsb_bern_delta_t(T& deg, T& sz, T& tot) {}
  T& deg_delta;
  T& sz_delta;
  T& tot_delta;
};

template <typename DT, typename DStat>
impl_increase_bottom_delta(DT& bottom,
			  const DStat& delta,
			  tag_bernoulli) {

  delta.sz_delta;
  delta.tot_delta;


  bottom.delta_stat_edge +=   delta.deg_delta;


  // stat.Nk
  // stat.Dik
  // stat.Zik
}

template <typename DT>
void impl_merge_left_right_delta(const DT& left, const DT& right, DT& intern,
                                 tag_bernoulli) {
  intern.delta_stat_edge = left.delta_stat_edge + right.delta_stat_edge;
  intern.delta_stat_total = left.delta_stat_total + right.delta_stat_total;
}

///////////////////////////////////////////
// Gam(dE + A) Gam(dH + B) Gam(A + B)	 //
// ----------------------- ------------- //
// Gam(dE + dH + A + B)    Gam(A) Gam(B) //
///////////////////////////////////////////
template <typename DT>
void impl_eval_delta_score(DT& D, tag_bernoulli) {
  static typename DT::Unit TOL = 1e-4;
  D.var_edge = D.stat_edge + D.a0 + TOL;
  D.var_hole = D.stat_total - D.stat_edge + D.b0 + TOL;

  // First term
  D.delta_score = fasterlgamma(D.delta_stat_edge + D.var_edge);
  D.delta_score +=
      fasterlgamma(D.delta_stat_total - D.delta_stat_edge + D.var_hole);
  D.delta_score -= fasterlgamma(D.delta_stat_total + D.var_edge + D.var_hole);

  // Second term
  D.delta_score += fasterlgamma(D.var_edge + D.var_hole);
  D.delta_score -= fasterlgamma(D.var_edge);
  D.delta_score -= fasterlgamma(D.var_hole);
}

///////////////////
// Gam(A + B)	 //
// ------------- //
// Gam(A) Gam(B) //
///////////////////
template <typename DT>
void impl_eval_score(DT& D, tag_bernoulli) {
  static typename DT::Unit TOL = 1e-4;
  D.var_edge = D.stat_edge + D.a0 + TOL;
  D.var_hole = D.stat_total - D.stat_edge + D.b0 + TOL;

  D.score = fasterlgamma(D.var_edge);
  D.score += fasterlgamma(D.var_hole);
  D.score -= fasterlgamma(D.var_edge + D.var_hole);
}

#endif
