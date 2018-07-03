#include "hsb_update_data.hh"

#ifndef HSB_FUNC_HH_
#define HSB_FUNC_HH_

template <typename Udata>
void remove_vertex(Udata& udata, const Index vertex_ii);

template <typename Udata>
void assign_vertex(Udata& udata, const Index vertex_ii, const Index cluster_kk);

////////////////////////////////////////////////////////////////

template <typename Udata>
void clear_tree_data(Udata& udata);

template <typename Udata>
void init_tree_var(Udata& udata);

template <typename Udata>
const auto& eval_tree_delta_score(Udata& udata, const Index ii);

template <typename Udata, typename Node>
const auto& eval_tree_delta_score(Udata& udata, Node r, const Index ii);

////////////////////////////////////////////////////////////////

template <typename Udata>
Scalar eval_tree_score(Udata& udata);

template <typename Udata>
Scalar eval_tree_stat_score(Udata& udata);

template <typename Udata>
void dump_tree_data_udata(Udata& udata);

template <typename Udata, typename Scalar>
void update_tree_var(Udata& udata, const Scalar rate);

////////////////////////////////////////////////////////////////

template <typename Udata, typename Node>
void increase_tree_stat(Udata& udata, Node r, const Index ii);

template <typename Udata>
void increase_tree_stat(Udata& udata);

template <typename Udata>
void increase_tree_stat(Udata& udata, const Index ii);

template <typename Udata, typename Node>
void decrease_tree_stat(Udata& udata, Node r, const Index ii);

template <typename Udata>
void decrease_tree_stat(Udata& udata);

template <typename Udata>
void decrease_tree_stat(Udata& udata, const Index ii);

#include "hsb_func_impl.hh"

#endif
