#ifndef _task_hpp_
#define _task_hpp_

#include <vector>

using namespace std;

class task;

class task {
 public:
  int ntask;
  vector<int> source_l_;
  vector<int> sink_l_;
  vector<int> delay_l_;
  vector<int> hop_l_;
  vector<int> cost_l_;

 public:
  task() {
    ;
  }
  task(int ntask) :
      ntask(ntask) {
    ;
  }
  int read_query(char *file_name, int graph_type) {
    char temp[255];
    FILE *topo_f;
    int id, src, snk, delay, hop, cost;
    topo_f = fopen(file_name, "r");
    vector<tuple<int, int, int, int>> task_l;
    if (graph_type == 1) {
      fgets(temp, 255, topo_f);
      while (fscanf(topo_f, "%d", &id) != EOF) {
        fscanf(topo_f, ",%d", &src);
        fscanf(topo_f, ",%d", &snk);
        fscanf(topo_f, ",%d", &delay);
        fscanf(topo_f, ",%d", &hop);
        fscanf(topo_f, ",%d", &cost);
        source_l_.push_back(src);
        sink_l_.push_back(snk);
        delay_l_.push_back(delay);
        hop_l_.push_back(hop);
        cost_l_.push_back(cost);
      }
    } else {
      int u, v, cost_lower, cost_upper, cost_limit;
      for (; fscanf(topo_f, "%d %d %d %d %d\n", &u, &v, &cost_lower, &cost_upper, &cost_limit) == 5;) {
        source_l_.push_back(u);
        sink_l_.push_back(v);
        delay_l_.push_back(cost_limit);
      }
    }
    fclose(topo_f);
    ntask = source_l_.size();

    return 0;
  }

  // used for COLA data set
  int load_query(char *file_name) {
    FILE *topo_f;
    topo_f = fopen(file_name, "r");
    vector<tuple<int, int, int, int>> task_l;
    int u, v, cost_lower, cost_upper, cost_limit;
    for (; fscanf(topo_f, "%d %d %d %d %d\n", &u, &v, &cost_lower, &cost_upper, &cost_limit) == 5;) {
      source_l_.push_back(u);
      sink_l_.push_back(v);
      delay_l_.push_back(cost_limit);
    }
    fclose(topo_f);
    ntask = source_l_.size();

    return 0;
  }
};

#endif
