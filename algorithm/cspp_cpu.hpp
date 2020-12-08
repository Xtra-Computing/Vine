#ifndef _cspp_cpu_hpp_
#define _cspp_cpu_hpp_
#include <vector>
#include <map>
#include "../include/label.hpp"
#include "../include/constraint.hpp"
#include "../include/edge.hpp"
#include "../include/graph.hpp"
#include "../include/util.hpp"
#include "../include/TimeMeasurer.hpp"
#include "../include/graph.hpp"
#include <omp.h>

#define LRG_FLOAT 1000000

template<typename Label_t>
void add_new_frontier(std::map<int, std::vector<Label_t>> &next_round, Label_t &new_label) {
  if (next_round.find(new_label.vid) != next_round.end()) {
    bool dominate = true;
    for (auto it = next_round[new_label.vid].begin(); it != next_round[new_label.vid].end();) {
      if (Label_t::dominate_cpu(*it, new_label)) {
        // new label is dominated
        dominate = false;
        break;
      }
      if (Label_t::dominate_cpu(new_label, *it)) {
        // new label dominates some other labels
        it = next_round[new_label.vid].erase(it);
      } else {
        it++;
      }
    }
    if (!dominate) {
      return; // go to next edge
    }
    next_round[new_label.vid].push_back(new_label);
  } else {
    std::vector<Label_t> new_frontier;
    new_frontier.push_back(new_label);
    next_round[new_label.vid] = new_frontier;
  }
}

template<typename Label_t, typename Edge_t>
int RUN_CSPP_CPU(
    Graph<Edge_t> *g,
    int ntask,
    vector<Label_t> &init_state_l,
    vector<int> &sink_l,
    vector<int *> &path_l,
    vector<int> &size_of_result_l,
    vector<c_constraint<Label_t>> &C_l) {

//====================
// timer starts
//====================

  TimeMeasurer total_time;
  total_time.tic();
  for (int i = 0; i < ntask; i++) {
//    std::cout << C_l[i].delay_L << std::endl;
    std::cout << "==================" << std::endl;
    // prepare data
    auto src = init_state_l[i];
    std::vector<Label_t> frontiers;
    std::map<int, std::vector<Label_t>> dominace_list;
    std::map<int, std::vector<Label_t>> next_round;
    frontiers.push_back(src);
    TimeMeasurer t;
    t.tic();
    {
      while (!frontiers.empty()) {
//      std::cout << "==================" << std::endl;
        for (int ff = 0; ff < frontiers.size(); ff++) {
          auto frontier = frontiers[ff];
          // expand each of it
          for (auto e : g->v_l[frontier.vid].edges) {
            Label_t new_label;
            frontier.expand_cpu(e->snk, *e, &new_label);
            if (!C_l[i].pass_constraints_check_cpu(new_label)) {
              continue;
            }
            if (dominace_list.find(new_label.vid) != dominace_list.end()) {
              bool dominate = true;
              for (auto it = dominace_list[new_label.vid].begin(); it != dominace_list[new_label.vid].end();) {
                if (Label_t::dominate_cpu(*it, new_label)) {
                  // new label is dominated
                  dominate = false;
                  break;
                }
                if (Label_t::dominate_cpu(new_label, *it)) {
                  // new label dominates some other labels
                  it = dominace_list[new_label.vid].erase(it);
                } else {
                  it++;
                }
              }
              if (!dominate) {
                continue; // go to next edge
              }
              dominace_list[new_label.vid].push_back(new_label);
              add_new_frontier(next_round, new_label);
            } else {
              std::vector<Label_t> new_dominace_vec;
              new_dominace_vec.push_back(new_label);
              dominace_list[new_label.vid] = new_dominace_vec;
              add_new_frontier(next_round, new_label);
            }
          }
        }
        frontiers.clear();
        for (auto i : next_round) {
          frontiers.insert(frontiers.end(), i.second.begin(), i.second.end());
        }
        next_round.clear();
      }
    }
    t.toc();
    t.print_ms("==");
  }
  total_time.toc();
  total_time.print_ms("Total Time is");

  return 0;
}

template<typename Label_t, typename Edge_t>
int RUN_CSPP_OMP(
    Graph<Edge_t> *g,
    int ntask,
    vector<Label_t> &init_state_l,
    vector<int> &sink_l,
    vector<int *> &path_l,
    vector<int> &size_of_result_l,
    vector<c_constraint<Label_t>> &C_l) {

//====================
// timer starts
//====================

  TimeMeasurer total_time;
  total_time.tic();
  for (int i = 0; i < ntask; i++) {
//    std::cout << C_l[i].delay_L << std::endl;
    std::cout << "==================" << std::endl;
    // prepare data
    auto src = init_state_l[i];
    std::vector<Label_t> frontiers;
    std::vector<Label_t> next_round;

    std::vector<std::vector<Label_t>> dominace_list;
    for (int j = 0; j < g->N; j++) {
      dominace_list.emplace_back(std::vector<Label_t>());
    }
    frontiers.push_back(src);
    TimeMeasurer t;
    t.tic();

    while (!frontiers.empty()) {
#pragma omp parallel for
      for (int ff = 0; ff < frontiers.size(); ff++) {
        auto frontier = frontiers[ff];
        // expand each of it
        for (auto e : g->v_l[frontier.vid].edges) {
          Label_t new_label;
          frontier.expand_cpu(e->snk, *e, &new_label);
          if (!C_l[i].pass_constraints_check_cpu(new_label)) {
            continue;
          }
          bool dominate = true;
          for (auto it = dominace_list[new_label.vid].begin(); it != dominace_list[new_label.vid].end();) {
            if (Label_t::dominate_cpu(*it, new_label)) {
              // new label is dominated
              dominate = false;
              break;
            }
            if (Label_t::dominate_cpu(new_label, *it)) {
              // new label dominates some other labels
#pragma omp critical
              it = dominace_list[new_label.vid].erase(it);
            } else {
              it++;
            }
          }
          if (!dominate) {
            continue; // go to next edge
          }
#pragma omp critical
          {
            dominace_list[new_label.vid].push_back(new_label);
            next_round.push_back(new_label);
          }
        }
      }
#pragma omp single
      {
        frontiers.swap(next_round);
        next_round.clear();
      }
#pragma omp barrier
    }

    for (auto i : dominace_list[g->N - 1]) {
      std::cout << i.d << " " << i.delay << std::endl;
    }
    t.toc();
    t.print_ms("==");
  }
  total_time.toc();
  total_time.print_ms("Total Time is");

  return 0;
}
#endif
