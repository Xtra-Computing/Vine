#ifndef _cspp_hpp_
#define _cspp_hpp_
#include <cuda_runtime.h>
#include <vector>
#include "../include/label.hpp"
#include "../include/constraint.hpp"
#include "../include/graph.hpp"
#include "cspp.cuh"
#include "cspp_cpu.hpp"

template<typename Label_t, typename Edge_t>
class cspp {
  using Graph_t = Graph<Edge_t>;
  using Constraint_t = c_constraint<Label_t>;
 public:
  int *path;

  vector<Label_t> init_state;
  vector<int> sink_l;
  int ntask;
  vector<int *> path_l;
  vector<int> size_of_result_l;

  Graph_t *g;
  vector<Constraint_t> &C_l;
 public:
  // one constraints
  cspp(Graph_t *g, Label_t source_state, int sink, Constraint_t C) : g(g) {
    ntask = 1;
    init_state.push_back(source_state);
    C_l.push_back(C);
    path = (int *) malloc(g->N * sizeof(int));
    path_l.push_back(path);
  }

  cspp(Graph_t *g, int ntask, vector<Label_t> &source_state, vector<int> &sink_l, vector<Constraint_t> &C_l) :
      g(g), ntask(ntask), init_state(source_state), sink_l(sink_l), C_l(C_l) {
    // allocate double size for the path
    long N = g->N;
    int *p = (int *) malloc(sizeof(int) * N * (long) ntask);
    path_l.reserve(ntask);
    for (long i = 0; i < ntask; i++) {
      path_l.push_back(&p[N * i]);
    }
    size_of_result_l.reserve(ntask);
  }

  int solve() {
    std::cout << "number of tasks " << ntask << std::endl;
    RUN_CSPP<Label_t, Edge_t>(g, ntask, init_state, sink_l, path_l, size_of_result_l, C_l);

    return 0;
  }

  int solve_cpu() {
    RUN_CSPP_CPU<Label_t, Edge_t>(g, ntask, init_state, sink_l, path_l, size_of_result_l, C_l);

    return 0;
  }

  int solve_omp() {
    RUN_CSPP_OMP<Label_t, Edge_t>(g, ntask, init_state, sink_l, path_l, size_of_result_l, C_l);
  }

  int self_validate() {

    return 0;
  }
};

#endif
