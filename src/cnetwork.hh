// (c) Yongjin Park, 2013
// collapsible network data structure
#include "rcpp_util.hh"
#include "sparse_data.hh"

#ifndef COLLAPSABLE_NETWORK_HH
#define COLLAPSABLE_NETWORK_HH

struct collapsible_network_t {
  // copy constructor
  explicit collapsible_network_t(const collapsible_network_t& C)
      : edges(C.edges),
        within_edges(C.within_edges),
        within_degrees(C.within_degrees),
        degree_sum(C.degree_sum),
        children(C.children),
        size_vec(C.size_vec),
        Dtot(C.Dtot) {}

  // just pre-calculate for faster access
  explicit collapsible_network_t(const sp_stat_mat_t& G0) : edges(G0) {
    // no within-edges, within-degrees, no children
    // degrees and total degree sum
    for (sp_stat_mat_t::row_iterator_t ri = G0.begin_row(); ri != G0.end_row();
         ++ri) {
      double d = sum(G0.pull_const_row(dim_t(*ri)), 1);
      degree_sum(dim_t(*ri), d);
      size_vec(dim_t(*ri), 1.);
    }
    Dtot = sum(degree_sum, 1);
  }

  // contract two nodes: u + v -> u
  void contract_pair(int u, int v);

 private:
  // delete node and associate data
  // to release non-existing node
  void clear_node(int u);

  // delete edge
  void clear_edge(int u, int v);

 public:
  sp_stat_mat_t edges;           // E[u,v] = edge weights (u,v)
  sp_stat_vec_t within_edges;    // E[u] = sum_{i<j} e_ij I[i,j in C[u]]
  sp_stat_vec_t within_degrees;  // wd[u] = sum_{i<j} d_i d_j I[i,j in C[u]]
  sp_stat_vec_t degree_sum;      // d[u] = degree sum of collapsed nodes
  sp_stat_mat_t children;        // C[u] = collapsed nodes under u
  sp_stat_vec_t size_vec;        // size = |C[u]| + 1
  double Dtot;                   // total sum of degrees

};  // End of collapsed network type

// collapse network G0 and return smaller network G
std::shared_ptr<collapsible_network_t> collapse_network(
    const collapsible_network_t& G0);

// matching node pairs by Jaccard coefficient
void jaccard_matching(const sp_stat_mat_t& G0, const sp_stat_vec_t& degree,
                      std::unordered_map<int, int>& Match);

////////////////////////////////////////////////////////////////
// error-checking after (many) collapsing / pair contraction
// comparing with the very original network
// - inefficient implementation
void check_collapsing_error(const collapsible_network_t& G0,
                            collapsible_network_t& G);

#endif  // End of collpased_network_hh
