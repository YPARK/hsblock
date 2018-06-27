#ifndef OPTIONS_HH_
#define OPTIONS_HH_

struct options_t {
  explicit options_t() {
    VBITER = 2000;
    K = 1;
    TREE_DEPTH = 2;
    VERBOSE = true;
    RSEED = 19937;
  }

  const int vbiter() const { return VBITER; };
  const int tree_depth() const { return TREE_DEPTH; };
  const int k() const { return K; };
  const int rseed() const { return RSEED; };
  const bool verbose() const { return VERBOSE; }

  int VBITER;
  int K;
  int TREE_DEPTH;
  int RSEED;

  bool VERBOSE;
};

#endif
