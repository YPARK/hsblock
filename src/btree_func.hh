#ifndef BTREE_FUNC_HH_
#define BTREE_FUNC_HH_

#include "btree.hh"

////////////////////////////////////////////////////////////////
// Typical statistical inference
// 1. Collect sufficient statistics (bottom -> top)
// 2. Calculate partial scores      (bottom -> top)
// 3. Sum up partial scores         (top -> down)


// A function for 

// E.g., a regular stochastic block model
//  Deg(G, i, r)
//  Tot(G, i, r)





// updater, propagator, 
// multiple b trees?

// just like score memoization


// stat[r->hash()] <= F.leaf(G, r)

// stat[r->hash()] <= F.lr(stat[r->left->hash()], stat[r->right->hash()])

// this stat can be Eigen matrix

// I can have multiple G'sx



template<typename D>
double max_score_path(typename btree_t<D>::node_ptr_t r,
                      double &deg_sum,
                      double &sz_sum,
                      int &argmax,
                      optim_state_t &ostate) {

  D &F = r->data;

  if (r->is_leaf()) {
    dim_t ii(ostate.ii);
    double d_i = ostate.G.degree_sum(ii);
    double Dtot = ostate.G.Dtot;

    int k = r->leaf_idx();
    dim_t kk(r->leaf_idx());
    deg_sum = ostate.clust_deg_mat(ii, kk);

#ifdef DEGREE_CORRECT
    sz_sum = d_i * ostate.vol(kk) / Dtot;
#else
    sz_sum = ostate.Bsz(kk);
#endif

#ifdef LCVI
    // // calculate log-predictive
    // double ret = F.log_pred(network_distrib_t::edge_t(deg_sum),
    //                         network_distrib_t::tot_t(sz_sum));
#endif
    argmax = k;
    return ret;
  }

  int arg_left, arg_right;
  double deg2left, deg2right, n2left, n2right;

  double Gleft =
      max_score_path(r->get_left(), deg2left, n2left, arg_left, ostate);
  double Gright =
      max_score_path(r->get_right(), deg2right, n2right, arg_right, ostate);

#ifdef LCVI
  // Gleft += F.log_pred(network_distrib_t::edge_t(deg2right),
  //                     network_distrib_t::tot_t(n2right));

  // Gright += F.log_pred(network_distrib_t::edge_t(deg2left),
  //                      network_distrib_t::tot_t(n2left));
#endif

  deg_sum = deg2left + deg2right;
  sz_sum = n2left + n2right;

  argmax = (Gleft > Gright) ? arg_left : arg_right;
  return (Gleft > Gright) ? Gleft : Gright;
}

////////////////////////////////////////////////////////////////

void calc_partial_score(model_tree_t::node_ptr_t r, double &deg_sum,
                        double &sz_sum, optim_state_t &ostate) {
  network_distrib_t &F = r->data;

  if (r->is_leaf()) {
    dim_t ii(ostate.ii);
    double d_i = ostate.G.degree_sum(ii);
    double Dtot = ostate.G.Dtot;

    dim_t k(r->leaf_idx());
    deg_sum = ostate.clust_deg_mat(ii, k);

#ifdef DEGREE_CORRECT
    sz_sum = d_i * ostate.vol(k) / Dtot;
#else
    sz_sum = ostate.Bsz(k);
#endif
    // cerr << deg_sum << ", " << sz_sum << ", " << w_e << ", " << w_d << endl;

#ifdef LCVI
    double score = F.log_pred(network_distrib_t::edge_t(deg_sum),
                              network_distrib_t::tot_t(sz_sum));
#else
    // calculate partial gradient
    double score = F.gradient(network_distrib_t::edge_t(deg_sum),
                              network_distrib_t::tot_t(sz_sum));
#endif
    ostate.left_score_memo(dim_t(r->hash()), score);
    return;
  }

  double deg2left, deg2right, n2left, n2right;

  calc_partial_score(r->get_left(), deg2left, n2left, ostate);
  calc_partial_score(r->get_right(), deg2right, n2right, ostate);

#ifdef LCVI
  double Gleft = F.log_pred(network_distrib_t::edge_t(deg2right),
                            network_distrib_t::tot_t(n2right));

  double Gright = F.log_pred(network_distrib_t::edge_t(deg2left),
                             network_distrib_t::tot_t(n2left));
#else
  double Gleft = F.gradient(network_distrib_t::edge_t(deg2right),
                            network_distrib_t::tot_t(n2right));

  double Gright = F.gradient(network_distrib_t::edge_t(deg2left),
                             network_distrib_t::tot_t(n2left));
#endif

  deg_sum = deg2left + deg2right;
  sz_sum = n2left + n2right;

  ostate.left_score_memo(dim_t(r->hash()), Gleft);
  ostate.right_score_memo(dim_t(r->hash()), Gright);
  return;
}

void sum_partial_score(model_tree_t::node_ptr_t r, sp_stat_vec_t &scores,
                       optim_state_t &ostate, double accum) {
  if (r->is_leaf()) {
    int k = r->leaf_idx();
    double memo = ostate.left_score_memo(dim_t(r->hash()));
    scores.increase(dim_t(k), accum + memo);
    return;
  }
  sum_partial_score(r->get_left(), scores, ostate,
                    accum + ostate.left_score_memo(dim_t(r->hash())));
  sum_partial_score(r->get_right(), scores, ostate,
                    accum + ostate.right_score_memo(dim_t(r->hash())));
}




#endif
