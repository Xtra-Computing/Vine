#ifndef _io_hpp_
#define _io_hpp_
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include "util.hpp"
#include "graph.hpp"
#include "task.hpp"
using namespace std;

extern float *A, *B, *Theta;

/**
 * @brief initialize graph
 * @tparam Graph_t
 * @param file_name
 * @param g
 * @param graph_type 1 for our graph, 2 for xxk graph
 * @return
 */
template<class Graph_t>
int read_graph(char *file_name, Graph_t *g, int graph_type) {
  char temp[255];
  FILE *topo_f;
  int id, src, snk;
  int weight, delay;
  topo_f = fopen(file_name, "r");
  int num_of_vertex, num_of_edge;
  if (graph_type == 1) {
    fscanf(topo_f, "%d", &num_of_vertex);
    fscanf(topo_f, "%c", temp);
    fscanf(topo_f, "%d", &num_of_edge);
    fscanf(topo_f, "%c", temp);
    num_of_edge *= 2;
    g->v_l.clear();
    g->v_l.reserve(num_of_vertex);
    g->N = num_of_vertex;
    g->e_l.clear();
    g->e_l.resize(num_of_edge);
    g->M = num_of_edge;
    using myVertex = Vertex<Graph_t>;
    for (int i = 0; i < num_of_vertex; i++) {
      myVertex *v = new myVertex(i);
      g->v_l.push_back(*v);
    }
    for (int i = 0; i < num_of_edge; i++) {
      g->e_l[i] = new Edge(i);
    }
    while (fscanf(topo_f, "%d", &src) != EOF) {
      fscanf(topo_f, "%d %d %d\n", &snk, &weight, &delay);
//      src--;
//      snk--;
      (g->e_l[id])->assign(id, src, snk, weight, delay);
      g->v_l[src].assign(g->e_l[id], g);
      id += 1;
      (g->e_l[id])->assign(id, snk, src, weight, delay);
      g->v_l[snk].assign(g->e_l[id], g);
      id += 1;
    }
  } else {
    fscanf(topo_f, "%d", &num_of_vertex);
    g->v_l.clear();
    g->v_l.reserve(num_of_vertex);
    g->N = num_of_vertex;
    g->e_l.clear();
    using myVertex = Vertex<Graph_t>;
    for (int i = 0; i < num_of_vertex; i++) {
      myVertex *v = new myVertex(i);
      g->v_l.push_back(*v);
    }
    int w, c, u, v, M = 0;
    for (; fscanf(topo_f, "%d:", &u) == 1;) {
      for (; fscanf(topo_f, "%d", &v) == 1;) {
        if (v == -1) break;
        int ret_val = fscanf(topo_f, "%d %d", &w, &c);
        if (ret_val < 0) {
          cout << "input file error" << endl;
          exit(0);
        }
        g->e_l.push_back(new Edge(M));
        g->e_l.back()->assign(M, u, v, w, c);
        g->v_l[u].assign(g->e_l[M++], g);
//        g->e_l.push_back(new Edge_Index(M));
//        g->e_l.back()->assign(M, v, u, w, c);
//        g->v_l[v].assign(g->e_l[M++], g);
      }
    }
    g->M = g->e_l.size();
    fclose(topo_f);
  }
  return 0;
}

#endif
