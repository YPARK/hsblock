#ifndef OPTIONS_HH_
#define OPTIONS_HH_

struct options_t {
  explicit options_t() {
    VBITER = 1000;
    INNER_ITER = 100;
    FINAL_INNER_ITER = 1000;
    BURNIN_ITER = 100;
    RECORD_INTERV = 10;
    K = 1;
    TREE_DEPTH = 2;
    VERBOSE = true;
    RSEED = 19937;
    DECAY = -0.55;
    DELAY = 1.0;
    RATE0 = 0.01;
  }

  const int vbiter() const { return VBITER; };
  const int burnin_iter() const { return BURNIN_ITER; };
  const int inner_iter() const { return INNER_ITER; };
  const int final_inner_iter() const { return FINAL_INNER_ITER; };
  const int record_interval() const { return RECORD_INTERV; };
  const int tree_depth() const { return TREE_DEPTH; };
  const int k() const { return K; };
  const int rseed() const { return RSEED; };
  const bool verbose() const { return VERBOSE; }
  const float decay() const { return DECAY; }
  const float delay() const { return DELAY; }
  const float rate0() const { return RATE0; }

  void set_inner_iter(const int new_iter) { INNER_ITER = new_iter; }

  int VBITER;
  int INNER_ITER;
  int FINAL_INNER_ITER;
  int BURNIN_ITER;
  int RECORD_INTERV;
  int K;
  int TREE_DEPTH;
  int RSEED;

  bool VERBOSE;
  float DECAY;
  float DELAY;
  float RATE0;
};

#endif
