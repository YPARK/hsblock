#include "cnetwork.hh"

void collapsible_network_t::contract_pair(int u, int v) {
  if (u == v) return;

#ifdef DEBUG
  assert(edges.has_row(dim_t(u)) && edges.has_row(dim_t(v)));
#endif

  // back-up (u,v) data and delete edge (u,v)
  double edge_uv = edges(dim_t(u), dim_t(v));
  double d_u = degree_sum(dim_t(u)), d_v = degree_sum(dim_t(v));

  // accumulate node data
  within_edges.increase(dim_t(u), edge_uv + within_edges(dim_t(v)));
  within_degrees.increase(dim_t(u), d_u * d_v + within_degrees(dim_t(v)));
  degree_sum.increase(dim_t(u), d_v);
  size_vec.increase(dim_t(u), size_vec(dim_t(v)));

  clear_edge(u, v);

  // accumulate edge data (v,w) + (u,w) -> (u,w)
  sp_stat_vec_t& adj_u = edges.pull_row(dim_t(u));
  sp_stat_vec_t& adj_v = edges.pull_row(dim_t(v));
  adj_u += adj_v;
  // copy row to column (ensure column version of adj_u correct)
  for (sp_stat_vec_t::iterator_t it = adj_u.begin_nonzero();
       it != adj_u.end_nonzero(); ++it) {
    int w = it.dim();

#ifdef DEBUG
    assert(w != v);  // can arise if there is self-loop
#endif

    edges(dim_t(w), dim_t(u), it.value());
  }

  // keep track of what are collapsed
  if (!children.has_row(dim_t(u))) children.add_row(dim_t(u));
  sp_stat_vec_t& under_u = children.pull_row(dim_t(u));
  // absorb children of v
  if (children.has_row(dim_t(v))) {
    under_u += children.pull_row(dim_t(v));
  }
  under_u(dim_t(v), 1);  // v is under u

  // delete node v in the end
  clear_node(v);
}

void collapsible_network_t::clear_node(int u) {
  // clear out adjacency
  const sp_stat_vec_t& adj_u = edges.pull_const_row(dim_t(u));
  for (sp_stat_vec_t::iterator_t it = adj_u.begin_nonzero();
       it != adj_u.end_nonzero(); ++it)
    if (it.dim() != u) edges(dim_t(it.dim()), dim_t(u), 0.);

  edges.del_row(dim_t(u));

  // clear out node data
  degree_sum(dim_t(u), 0.);
  within_edges(dim_t(u), 0.);
  within_degrees(dim_t(u), 0.);
  size_vec(dim_t(u), 0.);
  children.del_row(dim_t(u));
}

void collapsible_network_t::clear_edge(int u, int v) {
  edges(dim_t(u), dim_t(v), 0.);
  edges(dim_t(v), dim_t(u), 0.);
}

////////////////////////////////////////////////////////////////
// error-checking after (many) collapsing / pair contraction
// comparing with the very original network
// - inefficient implementation
void check_collapsing_error(const collapsible_network_t& G0,
                            collapsible_network_t& G) {
  const sp_stat_mat_t& edges0 = G0.edges;
  const sp_stat_vec_t& deg0 = G0.degree_sum;
  sp_stat_mat_t& edges = G.edges;
  sp_stat_mat_t& children = G.children;
  sp_stat_vec_t& within_edges = G.within_edges;
  sp_stat_vec_t& size_vec = G.size_vec;
  sp_stat_vec_t& within_degrees = G.within_degrees;
  sp_stat_vec_t& degsum = G.degree_sum;

  vector<int> set_r;
  vector<int> set_c;
  int r, c;
  for (sp_stat_mat_t::row_iterator_t ri = edges.begin_row();
       ri != edges.end_row(); ++ri) {
    r = *ri;
    set_r.clear();
    set_r.push_back(r);

    // cout << r;
    if (children.has_row(dim_t(r))) {
      const sp_stat_vec_t& r_vec = children.pull_const_row(dim_t(r));
      for (sp_stat_vec_t::iterator_t it = r_vec.begin_nonzero();
           it != r_vec.end_nonzero(); ++it) {
        set_r.push_back(it.dim());
        // cout << " " << it.dim();
      }
    }
    // cout << endl;
    if (abs(set_r.size() - size_vec(dim_t(r))) > 1e-5)
      cout << " r=" << r << " size " << set_r.size() << " vs "
           << size_vec(dim_t(r)) << endl;
    assert(abs(set_r.size() - size_vec(dim_t(r))) < 1e-5);

    // within edges and degrees
    if (set_r.size() == 1) {
      assert(abs(within_edges(dim_t(r))) < 1e-5);
      assert(abs(within_degrees(dim_t(r))) < 1e-5);
    } else {
      // collapsed edges
      double val_org = 0.;
      for (int i = 0; i < set_r.size(); ++i)
        for (int j = (i + 1); j < set_r.size(); ++j)
          val_org += edges0(dim_t(set_r[i]), dim_t(set_r[j]));

      double val_debug = within_edges(dim_t(r));
      if (abs(val_org - val_debug) > 1e-5)
        cout << " within edge " << r << " : " << val_org << " vs " << val_debug
             << endl;
      assert(abs(val_debug - val_org) < 1e-5);

      // degree product and sum
      val_org = 0.;
      for (int i = 0; i < set_r.size(); ++i)
        for (int j = (i + 1); j < set_r.size(); ++j)
          val_org += deg0(dim_t(set_r[i])) * deg0(dim_t(set_r[j]));

      val_debug = within_degrees(dim_t(r));
      if (abs(val_org - val_debug) > 1e-5)
        cout << " within degree " << r << " : " << val_org << " vs "
             << val_debug << endl;
      assert(abs(val_debug - val_org) < 1e-5);

      val_org = 0.;
      for (int i = 0; i < set_r.size(); ++i) val_org += deg0(dim_t(set_r[i]));

      val_debug = degsum(dim_t(r));
      if (abs(val_org - val_debug) > 1e-5)
        cout << " degree sum " << r << " : " << val_org << " vs " << val_debug
             << endl;
      assert(abs(val_debug - val_org) < 1e-5);
    }

    // between edges
    const sp_stat_vec_t& row = edges.pull_const_row(dim_t(r));
    for (sp_stat_vec_t::iterator_t ci = row.begin_nonzero();
         ci != row.end_nonzero(); ++ci) {
      c = ci.dim();
      set_c.clear();
      set_c.push_back(c);
      // cout << c;
      if (children.has_row(dim_t(c))) {
        const sp_stat_vec_t& c_vec = children.pull_const_row(dim_t(c));
        for (sp_stat_vec_t::iterator_t it = c_vec.begin_nonzero();
             it != c_vec.end_nonzero(); ++it) {
          set_c.push_back(it.dim());
          // cout << " " << it.dim();
        }
      }
      // cout << endl;

      // calculate from original edges
      double val_org = 0.;
      for (int i = 0; i < set_r.size(); ++i)
        for (int j = 0; j < set_c.size(); ++j)
          val_org += edges0(dim_t(set_r[i]), dim_t(set_c[j]));

      double val_debug = edges(dim_t(r), dim_t(c));
      if (abs(val_org - val_debug) > 1e-5)
        cout << " edge " << r << ", " << c << " : " << val_org << " vs "
             << val_debug << endl;
      assert(abs(val_org - val_debug) <= 1e-5);
    }
  }
}
