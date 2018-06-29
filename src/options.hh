#ifndef OPTIONS_HH_
#define OPTIONS_HH_

struct options_t {
  explicit options_t() {
    VBITER = 2000;
    INNER_ITER = 1000;
    BURNIN_ITER = 100;
    RECORD_INTERV = 10;
    K = 1;
    TREE_DEPTH = 2;
    VERBOSE = true;
    ECONOMY = false;
    RSEED = 19937;
    DECAY = -0.55;
    RATE0 = 0.01;
  }

  const int vbiter() const { return VBITER; };
  const int burnin_iter() const { return BURNIN_ITER; };
  const int inner_iter() const { return INNER_ITER; };
  const int record_interval() const { return RECORD_INTERV; };
  const int tree_depth() const { return TREE_DEPTH; };
  const int k() const { return K; };
  const int rseed() const { return RSEED; };
  const bool verbose() const { return VERBOSE; }
  const bool economy() const { return ECONOMY; }
  const float decay() const { return DECAY; }
  const float delay() const { return DELAY; }
  const float rate0() const { return RATE0; }

  int VBITER;
  int INNER_ITER;
  int BURNIN_ITER;
  int RECORD_INTERV;
  int K;
  int TREE_DEPTH;
  int RSEED;

  bool VERBOSE;
  bool ECONOMY;

  float DECAY;
  float DELAY;
  float RATE0;
};

#endif
