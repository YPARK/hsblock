////////////////////////////////////////////////////////////////
// a simple binary tree class
// (c) Yongjin Park, 2012

#include <cmath>
#include <memory>
#include <unordered_map>
#include <vector>
#include "rcpp_util.hh"

#ifndef BTREE_HH_
#define BTREE_HH_

////////////////////////////////////////////////////////////////
// To find trailing zeros efficiently
// (source http://graphics.stanford.edu/~seander/bithacks.html)
static const int MultiplyDeBruijnBitPosition[32] = {
    0,  1,  28, 2,  29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4,  8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6,  11, 5,  10, 9};

#define tree_return_false(msg)                           \
  {                                                      \
    Rcpp::Rcerr << "[TREE] Error: " << msg << std::endl; \
    return false;                                        \
  }
#define tree_throw(msg)                                  \
  {                                                      \
    Rcpp::Rcerr << "[TREE] Error: " << msg << std::endl; \
    exit(1);                                             \
  }
#define tree_message(msg) \
  { Rcpp::Rcerr << "[TREE] " << msg << std::endl; }

////////////////////////////////////////////////////////////////
// fixed binary tree class
template <typename D>
class btree_t {
 public:
  typedef D node_data_t;
  typedef std::shared_ptr<node_data_t> node_data_ptr_t;

 private:
  // we use conventional binary heap-like tree
  // constant time easy access
  // internal index for tree nodes

  // **************** useful facts ****************
  // root index = 0 and a single root is depth=1
  // clearly depth > 1, otherwise useless
  // pa[i] = floor( (i-1)/2 )
  // l[i] = 2*i + 1
  // r[i] = 2*i + 2
  // ntot = 2^depth -1
  // nleaves = 2^depth
  // leaf index = [2^depth-1 .. 2^(1+depth)-1)
  // 0-based indexing

  int depth;        // depth of a tree
  int ntot;         // total number of tree nodes
  int leaf_offset;  // (1 << (depth-1)) - 1

  // For faster lookups of a LCA
  // we use in-order traversal index
  // LCA(x,y) could be done in O(1)

  std::vector<int> inorder2tree;  // inorder idx -> tree idx
  std::vector<int> tree2inorder;  // tree idx -> inorder idx
  std::vector<node_data_ptr_t> node_data_map;

  // private helper functions
  void init();
  void clear();

  // Recursively build a inorder traversal
  void dfs_inorder_idx(int root, int& idx);

  int left_idx(int ii) const;
  int right_idx(int ii) const;
  int pa_idx(int ii) const;
  int leaf_internal_idx(int k) const;

  // find the number of trailing zeros in 32-bit v
  int count_zeros_right(int v) const;

  // Fast access of lowest common ancestor of leaf_x and leaf_y
  // implemetation of Moufatich (2008)
  // takes internal index of (x, y) in [ leaf_offset .. ntot )
  // outputs internal index in [ 0 .. num_internals() )
  int get_lca(int xi, int yi) const;

  // lca node idx
  int get_lca_idx(int leaf_x, int leaf_y) const;

 public:
  explicit btree_t(int);
  virtual ~btree_t();

  // leaf_x, leaf_y in [0 .. #leaves)
  node_data_t& get_lca_node(int leaf_x, int leaf_y);

  // leaf node data
  node_data_t& get_leaf_node(int leaf);

  // n-th node
  node_data_t& get_nth_node(int n);

  // leaf_x, leaf_y in [0 .. #leaves)
  const node_data_t& get_lca_node_const(int leaf_x, int leaf_y) const;

  // leaf node data
  const node_data_t& get_leaf_node_const(int leaf) const;

  // n-th node
  const node_data_t& get_nth_node_const(int n) const;

  int num_leaves() const;
  int num_internals() const;
  int num_nodes() const;
  int get_depth() const;

 public:
  // a wrapper class for node idx and data
  struct node_t {
    int id;
    node_data_t& data;
    btree_t& tree;

    node_t(int _id, node_data_t& _nd, btree_t& _t) : data(_nd), tree(_t) {
      id = _id;
    }

    node_t(node_t& other) : data(other.data), tree(other.tree) {
      id = other.id;
    }

    int leaf_idx() { return (id - tree.leaf_offset); }

    bool is_leaf() {
      return (leaf_idx() >= 0 && leaf_idx() < tree.num_leaves());
    }
    bool has_pa() { return id > 0 && tree.pa_idx(id) >= 0; }

    const int hash() { return id; }

    std::shared_ptr<node_t> get_pa() {
#ifdef DEBUG
      assert(has_pa());
#endif

      int id_pa = tree.pa_idx(id);
      return tree.node_obj_map[id_pa];
    }

    std::shared_ptr<node_t> get_left() {
#ifdef DEBUG
      assert(!is_leaf());
#endif

      int id_left = tree.left_idx(id);
      return tree.node_obj_map[id_left];
    }

    std::shared_ptr<node_t> get_right() {
#ifdef DEBUG
      assert(!is_leaf());
#endif

      int id_right = tree.right_idx(id);
      return tree.node_obj_map[id_right];
    }

    friend class btree_t;
  };

 public:
  friend class node_t;

  typedef std::shared_ptr<node_t> node_ptr_t;

  node_ptr_t root_node_obj();
  node_ptr_t left_node_obj(node_ptr_t node);
  node_ptr_t right_node_obj(node_ptr_t node);
  node_ptr_t pa_node_obj(node_ptr_t node);

  // leaf_x, leaf_y in [0 .. #leaves)
  node_ptr_t get_lca_node_obj(int leaf_x, int leaf_y);

  // just n-th one
  node_ptr_t get_nth_node_obj(int n);

  // leaf node
  node_ptr_t get_leaf_node_obj(int leaf);

 private:
  // pre-generated map
  // node idx to node pointer
  std::vector<std::shared_ptr<node_t> > node_obj_map;
};

///////////////////////////////////////////////////////////////
template <typename D>
void btree_t<D>::init() {
  ntot = (1 << (1 + depth)) - 1;
  inorder2tree.resize(ntot, 0);
  tree2inorder.resize(ntot, 0);
  node_data_map.clear();
  leaf_offset = (1 << depth) - 1;
  int idx = 0;              // cumulative index
  dfs_inorder_idx(0, idx);  // build look-ups

  for (int j = 0; j < ntot; ++j) {
    node_data_ptr_t node_data_ptr(new node_data_t);
    node_ptr_t node_obj_ptr(new node_t(j, *(node_data_ptr.get()), *this));
    node_data_map.push_back(node_data_ptr);
    node_obj_map.push_back(node_obj_ptr);
  }
#ifdef DEBUG
  assert(ntot == node_data_map.size());
  assert(ntot == inorder2tree.size());
  assert(ntot == tree2inorder.size());
  tree_message("Initialized");
#endif

  return;
}

template <typename D>
void btree_t<D>::clear() {}

////////////////////////////////////////////////////////////////
// Recursively build a inorder traversal
template <typename D>
void btree_t<D>::dfs_inorder_idx(int root, int& idx) {
  int left = left_idx(root);
  int right = right_idx(root);
  if (std::min(left, right) >= ntot) {  // can't go further down
    inorder2tree[idx] = root;
    tree2inorder[root] = idx;
    ++idx;
    return;
  }
  dfs_inorder_idx(left, idx);
  inorder2tree[idx] = root;
  tree2inorder[root] = idx;
  ++idx;
  dfs_inorder_idx(right, idx);

  return;
}

template <typename D>
int btree_t<D>::left_idx(int ii) const {
  return 2 * ii + 1;
}

template <typename D>
int btree_t<D>::right_idx(int ii) const {
  return 2 * ii + 2;
}

template <typename D>
int btree_t<D>::pa_idx(int ii) const {
  return (ii - 1) / 2;
}

template <typename D>
int btree_t<D>::leaf_internal_idx(int k) const {
#ifdef DEBUG
  int ret = (k + leaf_offset);
  assert(k >= 0 && ret < ntot);
  return ret;
#else

  return (k + leaf_offset);
#endif
}

////////////////////////////////////////////////////////////////
// find the number of trailing zeros in 32-bit v
template <typename D>
int btree_t<D>::count_zeros_right(int v) const {
  int r;  // result goes here
  r = MultiplyDeBruijnBitPosition[((uint32_t)((v & -v) * 0x077CB531U)) >> 27];
  return r;
}

////////////////////////////////////////////////////////////////
// Fast access of lowest common ancestor of leaf_x and leaf_y
// implemetation of Moufatich (2008)
// takes internal index of (x, y) in [ leaf_offset .. ntot )
// outputs internal index in [ 0 .. num_internals() )
template <typename D>
int btree_t<D>::get_lca(int xi, int yi) const {
#ifdef DEBUG
  assert(xi != yi);
#endif

  int x, y, idx1, idx2, idx3, idx, lca;
  x = tree2inorder[xi] + 1;  // inorder index (1-based)
  y = tree2inorder[yi] + 1;  // inorder index (1-based)

  int xy_xor = x ^ y;           // take difference of two bits
  idx1 = (int)log2(xy_xor);     // position of leftmost 1-bit
  idx2 = count_zeros_right(x);  // #0's until rightmost 1-bit
  idx = std::max(idx1, idx2);   //
  idx3 = count_zeros_right(y);  // #0's until rightmost 1-bit
  if (idx3 > idx) {
    lca = y;
  } else {
    lca = (x >> (idx + 1)) << (idx + 1);
  }
  lca = (lca | (1 << idx));

  return inorder2tree.at(lca - 1);
}

////////////////////////////////////////////////////////////////
// lca node idx
template <typename D>
int btree_t<D>::get_lca_idx(int leaf_x, int leaf_y) const {
  if (leaf_x == leaf_y) return leaf_internal_idx(leaf_x);

  return get_lca(leaf_internal_idx(leaf_x), leaf_internal_idx(leaf_y));
}

template <typename D>
btree_t<D>::btree_t(int _depth) {
  depth = std::max(1, _depth);  // depth should be at least 2
  init();
}

template <typename D>
btree_t<D>::~btree_t() {
  clear();
}

////////////////////////////////////////////////////////////////
// leaf_x, leaf_y in [0 .. #leaves)
template <typename D>
typename btree_t<D>::node_data_t& btree_t<D>::get_lca_node(int leaf_x,
                                                           int leaf_y) {
  // int lca = get_lca( leaf_internal_idx(leaf_x), leaf_internal_idx(leaf_y) );
  // tree_message( "has lca: " << lca << " of x=" << leaf_internal_idx(leaf_x)
  // << " y=" << leaf_internal_idx(leaf_y) );
  int lca = get_lca_idx(leaf_x, leaf_y);
  return *(node_data_map[lca].get());
}

// leaf node data
template <typename D>
typename btree_t<D>::node_data_t& btree_t<D>::get_leaf_node(int leaf) {
  return *(node_data_map[leaf_internal_idx(leaf)].get());
}

//
template <typename D>
typename btree_t<D>::node_data_t& btree_t<D>::get_nth_node(int n) {
#ifdef DEBUG
  assert(n >= 0 && n < ((int)ntot));
#endif

  return *(node_data_map[n].get());
}

// leaf_x, leaf_y in [0 .. #leaves)
template <typename D>
const typename btree_t<D>::node_data_t& btree_t<D>::get_lca_node_const(
    int leaf_x, int leaf_y) const {
  // int lca = get_lca( leaf_internal_idx(leaf_x), leaf_internal_idx(leaf_y) );
  // tree_message( "has lca: " << lca << " of x=" << leaf_internal_idx(leaf_x)
  // << " y=" << leaf_internal_idx(leaf_y) );
  int lca = get_lca_idx(leaf_x, leaf_y);
  return *(node_data_map.at(lca).get());
}

// leaf node data
template <typename D>
const typename btree_t<D>::node_data_t& btree_t<D>::get_leaf_node_const(
    int leaf) const {
  return *(node_data_map.at(leaf_internal_idx(leaf)).get());
}

//
template <typename D>
const typename btree_t<D>::node_data_t& btree_t<D>::get_nth_node_const(
    int n) const {
#ifdef DEBUG
  assert(n >= 0 && n < ((int)ntot));
#endif

  return *(node_data_map.at(n).get());
}

template <typename D>
int btree_t<D>::num_leaves() const {
  return (ntot - leaf_offset);
}

template <typename D>
int btree_t<D>::num_internals() const {
  return leaf_offset;
}

template <typename D>
int btree_t<D>::num_nodes() const {
#ifdef DEBUG
  assert(ntot == node_data_map.size());
  assert(ntot == inorder2tree.size());
  assert(ntot == tree2inorder.size());
#endif

  return ntot;
}

template <typename D>
int btree_t<D>::get_depth() const {
  return depth;
}

////////////////////////////////////////////////////////////////
template <typename D>
typename btree_t<D>::node_ptr_t btree_t<D>::root_node_obj() {
  // return node_ptr_t( new node_t(0,*(node_data_map[0].get()),*this) );
  return node_obj_map[0];
}

template <typename D>
typename btree_t<D>::node_ptr_t btree_t<D>::left_node_obj(node_ptr_t node) {
#ifdef DEBUG
  assert(!(node->is_leaf()));
#endif

  int lid = left_idx(node->id);
  // return node_ptr_t( new node_t(lid,*(node_data_map[lid].get()),*this) );
  return node_obj_map[lid];
}

template <typename D>
typename btree_t<D>::node_ptr_t btree_t<D>::right_node_obj(node_ptr_t node) {
#ifdef DEBUG
  assert(!(node->is_leaf()));
#endif

  int rid = right_idx(node->id);
  //   return node_ptr_t( new node_t(rid,*(node_data_map[rid].get()),*this) );
  return node_obj_map[rid];
}

template <typename D>
typename btree_t<D>::node_ptr_t btree_t<D>::pa_node_obj(node_ptr_t node) {
#ifdef DEBUG
  assert(node->has_pa());
#endif

  int pa = pa_idx(node->id);
  //   return node_ptr_t( new node_t(pa,*(node_data_map[pa].get()),*this) );
  return node_obj_map[pa];
}

// leaf_x, leaf_y in [0 .. #leaves)
template <typename D>
typename btree_t<D>::node_ptr_t btree_t<D>::get_lca_node_obj(int leaf_x,
                                                             int leaf_y) {
  // int lca = get_lca( leaf_internal_idx(leaf_x), leaf_internal_idx(leaf_y) );
  // tree_message( "has lca: " << lca << " of x=" << leaf_internal_idx(leaf_x)
  // << " y=" << leaf_internal_idx(leaf_y) );
  int lca = get_lca_idx(leaf_x, leaf_y);
  return node_obj_map[lca];
}

template <typename D>
typename btree_t<D>::node_ptr_t btree_t<D>::get_nth_node_obj(int n) {
#ifdef DEBUG
  assert(n >= 0 && n < ((int)ntot));
#endif

  return node_obj_map[n];
}

// leaf node data
template <typename D>
typename btree_t<D>::node_ptr_t btree_t<D>::get_leaf_node_obj(int leaf) {
  return node_obj_map[leaf_internal_idx(leaf)];
}

#endif /* BTREE_HH_ */
