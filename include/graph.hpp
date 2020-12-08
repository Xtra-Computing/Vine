#ifndef _graph_hpp_
#define _graph_hpp_
#include <cstdio>
#include <vector>
#include <cassert>
#include <algorithm>
#include <queue>
#include <set>
#include <list>
#include <limits>
#include <iostream>
#include <string>
#include <map>

#include "gpu_graph.hpp"
#include "edge.hpp"

using namespace std;
template<typename Edge_t>
class Graph;

template<typename Edge_t>
class gpuGraph;

class Edge {
 public:
  int id;
  int src, snk;
  float len, delay;
 public:
  Edge() {
    ;
  }
  Edge(int id) :
      id(id) {
    ;
  }
  Edge(int gg, int v, int I);

  int assign(
      int id,
      int src,
      int snk,
      float len = 0.0,
      float delay = 0.0) {
    this->id = id;
    this->src = src;
    this->snk = snk;
    this->len = len;
    this->delay = delay;

    return 0;
  }

  int getOwner() const {
    return src;
  }
  int operator()() const;
  int getVersion() const;
  bool operator<(const Edge &o) const;
};

template<typename Graph_t>
class Vertex {
 public:
  int id;
  vector<Edge *> edges;
  vector<int> next;
  vector<Edge *> edges_r;
  vector<int> pre;
 public:
  Vertex() {
    ;
  }
  Vertex(int id) :
      id(id) {
    ;
  }
  int assign(Edge *e, Graph_t *g) {
    this->edges.push_back(e);
    this->next.push_back(e->snk);
    Vertex *son = &(g->v_l[e->snk]);
    son->edges_r.push_back(e);
    son->pre.push_back(this->id);

    return 0;
  }
};

template<typename Edge_t>
class Graph {
 public:

  using MyV = Vertex<Graph<Edge_t> >;
  vector<MyV> v_l;
  vector<Edge *> e_l;
  int N, M;

  // CSR part
  // set as public for simpler manipulatign
  int *_count_row; // CSR pointer
  int *_eid_col; // CSR edge id
  int *_vid_col; //CSR tail v id
  int *e_sink; // eid to sink
  int *e_csr;  // eid to csr pos
  float *_len; // CSR edge len
  float *_delay; // CSR edge delay

  Edge_t *el;

  gpuGraph<Edge_t> *gg; // gpu graph on device
  gpuGraph<Edge_t> *gg_h; // gpu graph on host

 public:
  Graph() {
    // data allocation is outside... since public
    ;
  }
  Graph(int N, int M) :
      N(N), M(M), e_l(M) {
    v_l.reserve(N);
    e_l.reserve(2 * M);
    for (int i = 0; i < N; i++) {
      MyV *v = new MyV(i);
      v_l.push_back(*v);
    }
    for (int i = 0; i < M; i++) {
      e_l[i] = new Edge(i);
    }
  }

  int adj_to_csr();
  int make_gpu_graph();
};

template<typename Edge_t>
int Graph<Edge_t>::adj_to_csr() {
  // allocate memory for CSR
  _count_row = (int *) malloc(sizeof(int) * (N + 1));
  _eid_col = (int *) malloc(sizeof(int) * M);
  _vid_col = (int *) malloc(sizeof(int) * M);
  e_sink = (int *) malloc(sizeof(int) * M);
  e_csr = (int *) malloc(sizeof(int) * M);
  _len = (float *) malloc(sizeof(float) * M);
  _delay = (float *) malloc(sizeof(float) * M);

  el = (Edge_t *) malloc(sizeof(Edge_t) * M);
  _count_row[0] = 0;
  for (int i = 0; i < N; i++) {
    _count_row[i + 1] = _count_row[i] + v_l[i].edges.size();
    if (v_l[i].edges.empty()) {
      continue;
    }
    Edge *e;
    for (int j = _count_row[i], t = 0; j < _count_row[i + 1]; j++, t++) {
      e = v_l[i].edges[t];
      e_csr[e->id] = j;
      e_sink[e->id] = e->src;
      _eid_col[j] = e->id;
      _vid_col[j] = e->snk;
      _len[j] = e->len;
      _delay[j] = e->delay;
    }
  }

  return 0;
}

template<typename Edge_t>
int Graph<Edge_t>::make_gpu_graph() {
  gg_h = new gpuGraph<Edge_t>();
  gg_h->N = N;
  cudaMalloc(&(gg_h->_count_row_d), sizeof(int) * (N + 1));
  cudaMalloc(&(gg_h->_eid_col_d), sizeof(int) * M);
  cudaMalloc(&(gg_h->_vid_col_d), sizeof(int) * M);
  cudaMalloc(&(gg_h->_sink_d), sizeof(int) * M);
  cudaMalloc(&(gg_h->_csr_d), sizeof(int) * M);
  cudaMalloc(&(gg_h->_len_data_d), sizeof(float) * M);
  cudaMalloc(&(gg_h->_delay_data_d), sizeof(float) * M);
  cudaMalloc(&(gg_h->el_d), sizeof(Edge_t) * M);

  cudaMemcpy(gg_h->_count_row_d, _count_row, sizeof(int) * (N + 1), cudaMemcpyHostToDevice);
  cudaMemcpy(gg_h->_eid_col_d, _eid_col, sizeof(int) * M, cudaMemcpyHostToDevice);
  cudaMemcpy(gg_h->_vid_col_d, _vid_col, sizeof(int) * M, cudaMemcpyHostToDevice);
  cudaMemcpy(gg_h->_sink_d, e_sink, sizeof(int) * M, cudaMemcpyHostToDevice);
  cudaMemcpy(gg_h->_csr_d, e_csr, sizeof(int) * M, cudaMemcpyHostToDevice);
  cudaMemcpy(gg_h->_len_data_d, _len, sizeof(float) * M, cudaMemcpyHostToDevice);
  cudaMemcpy(gg_h->_delay_data_d, _delay, sizeof(float) * M, cudaMemcpyHostToDevice);
  cudaMemcpy(gg_h->el_d, el, sizeof(Edge_t) * M, cudaMemcpyHostToDevice);

  cudaMalloc(&(this->gg), sizeof(gpuGraph<Edge_t>));
  cudaMemcpy(this->gg, this->gg_h, sizeof(gpuGraph<Edge_t>), cudaMemcpyHostToDevice);

  return 0;
}

#endif
