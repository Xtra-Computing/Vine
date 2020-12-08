#include <stdio.h>
#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include "include/comm.hpp"
#include "include/task.hpp"
#include "include/label.hpp"
#include "include/constraint.hpp"
#include "include/label.hpp"

#ifndef CSP
#define CSP 1
#endif

using namespace std;

static Graph<c_edge> *g;

int pre_process_task(task *tasks, int N) { // to make sure source are within range
  for (int i = 0; i < tasks->ntask; i++) {
    tasks->source_l_[i] %= N;
    tasks->sink_l_[i] %= N;
  }

  return 0;
}

c_edge *create_edge_list(Graph<c_edge> *g) {
  c_edge *l = (c_edge *) malloc(sizeof(c_edge) * g->M);
  for (int i = 0; i < g->M; i++) {
    l[i].len = g->e_l[i]->len;
    l[i].delay = g->e_l[i]->delay;
  }

  return l;
}

int main(int argc, char *argv[]) {
  using Constraint = c_constraint<c_label>;
  using Graph = Graph<c_edge>;
  int graph_type, ntask;
  char graph_f[255], query_f[255];
  FILE *conf;
  conf = fopen("../inputs/config", "r");
  fscanf(conf, "%d", &graph_type);
  fscanf(conf, "%s", graph_f);
  fscanf(conf, "%s", query_f);
  fscanf(conf, "%d", &ntask);

  cudaSetDevice(0);
  g = new Graph();
  task *tasks = new task();
  read_graph(graph_f, g, graph_type);
  if (graph_type == 1) {
    tasks->read_query(query_f, 1);
    pre_process_task(tasks, g->N);
  } else {
    // only for COLA data set
    tasks->load_query(query_f);
  }
  g->adj_to_csr();
  g->el = create_edge_list(g);
  g->make_gpu_graph();

  cspp<c_label, c_edge> *solver;

  // create constraints, and source_l
  vector<Constraint> C_l;
  vector<c_label> source_state_l;
  for (int i = 0; i < ntask; i++) {
    c_label l(tasks->source_l_[i]);
    source_state_l.push_back(l);
    Constraint c(tasks->delay_l_[i]);
    C_l.push_back(c);
  }
  cuda_err_chk(cudaGetLastError());
  std::cout << "Number of V " << g->N << std::endl;
  std::cout << "Number of E " << g->M << std::endl;
  solver = new cspp<c_label, c_edge>(g, ntask, source_state_l, tasks->sink_l_, C_l);
  cuda_err_chk(cudaGetLastError());
  solver->solve();
  cuda_err_chk(cudaGetLastError());
  return 0;
}
