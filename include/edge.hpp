#ifndef _edge_hpp_
#define _edge_hpp_

template<typename Edge_t>
class Graph;

class edge {
 public:
  __device__ edge() {
  }

  ~edge() {
  }
};

class c_edge {
 public:
  float len;
  float delay;
 public:
  c_edge() {
  }
  c_edge(float l, float d) :
      len(l), delay(d) {
  }
};

#endif
