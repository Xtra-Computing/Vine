#ifndef _cspp_cuh_
#define _cspp_cuh_

#include <cuda_runtime.h>
#include "../include/label.hpp"
#include "../include/constraint.hpp"
#include "../include/edge.hpp"
#include "../include/graph.hpp"
#include "../include/util.hpp"
#include "thrust/scan.h"
#include "thrust/execution_policy.h"
#include "thrust/device_ptr.h"
#include "bucket.cuh"
#include "../include/TimeMeasurer.hpp"
#include "../include/AR.hpp"
#include <math.h>

#include <cuda_profiler_api.h>

#define USE_AR

#define LRG_FLOAT 1000000
#define NUM_SM 56
#define ITER 1000
#ifdef USE_AR
#define AR_DEGREE 3
#define THETA 0.95
#else
#define MAX_LABELS_IN_FRONTIER 1000000
#endif
int v_per_block = (BLOCK_SIZE / 32) * VSIZE;

template<typename Label_t>
double calculate_p(double ar_n, int n, int max_labels) {
  // M / size = N - P - theta*P + AR(n)P/n
  double est = max_labels * 1.0 / sizeof(Label_t) / (ar_n / n - THETA);
  if (est < 0) {
    return (double) max_labels;
  }
  return min(est, (double) max_labels);
}

int max_frontier(int totalLabel) {
  size_t free, total;
  cudaMemGetInfo(&free, &total);
  return free / sizeof(c_label) - totalLabel;
}

// 2 Cons massive with check function
template<typename Label_t, typename Edge_t>
int RUN_CSPP(
    Graph<Edge_t> *g,
    int ntask,
    vector<Label_t> &init_state_l,
    vector<int> &sink_l,
    vector<int *> &path_l,
    vector<int> &size_of_result_l,
    vector<c_constraint<Label_t>> &C_l);

template<typename Label_t>
__global__ void ck_CSPP_path(int source, Label_t **dominance_list, int *list_len, int *path, int *mutex);

template<typename Label_t>
__global__ void ck_init(int *dominance_list_len, Label_t **dominance_list, Label_t source) {
  dominance_list_len[source.vid] = 1;
  Label_t::copy(dominance_list[source.vid], &source);
}

template<typename Label_t, typename Edge_t>
int RUN_CSPP(
    Graph<Edge_t> *g,
    int ntask,
    vector<Label_t> &init_state_l,
    vector<int> &sink_l,
    vector<int *> &path_l,
    vector<int> &size_of_result_l,
    vector<c_constraint<Label_t>> &C_l) {
  int N = g->N;
  TimeMeasurer t;
  cuda_err_chk(cudaGetLastError());

  Label_t *buffer;
  cudaMallocManaged(&buffer, sizeof(Label_t) * N * SearchSpace);
  Label_t *bucket;
  cudaMallocManaged(&bucket, sizeof(Label_t) * N * SearchSpace);
  cuda_err_chk(cudaGetLastError());

  int *scan_out, *status;
  cudaMallocManaged(&scan_out, sizeof(int) * N * SearchSpace);
  cudaMallocManaged(&status, sizeof(int) * N * SearchSpace);

  int *mutex;
  cudaMallocManaged(&mutex, sizeof(int) * N * SearchSpace);
  cuda_err_chk(cudaGetLastError());

  cudaMemset(mutex, 0, sizeof(int) * N * SearchSpace);
  Label_t **dominance_list_d;
  cudaMallocManaged(&dominance_list_d, sizeof(Label_t *) * N);
  Label_t **dynamic_list_d;
  cudaMallocManaged(&dynamic_list_d, sizeof(Label_t *) * N / _K);
  cudaDeviceSynchronize();
  cuda_err_chk(cudaGetLastError());

  int *dominance_list_len;
  cudaMallocManaged(&dominance_list_len, sizeof(int) * N);
  int *dynamic_list_len;
  cudaMallocManaged(&dynamic_list_len, sizeof(int) * N / _K);

  Label_t *dominance_list;
  Label_t *dynamic_list;
  vector<Label_t *> dominance_list_l;
  vector<Label_t *> dynamic_list_l;
  cudaMallocManaged(&dominance_list, sizeof(Label_t) * N * SearchSpace);
  cudaMallocManaged(&dynamic_list, sizeof(Label_t) * N * SearchSpace * GlobalSpaceFactor);

  for (int T = 0; T < N; T++) {
    dominance_list_l.push_back(&dominance_list[T * SearchSpace]);
    if (T % _K == 0)
      dynamic_list_l.push_back(&dynamic_list[T * SearchSpace * GlobalSpaceFactor]);
  }
  cudaMemcpy(dominance_list_d, &dominance_list_l[0], sizeof(Label_t *) * N, cudaMemcpyHostToDevice);
  cudaMemcpy(dynamic_list_d, &dynamic_list_l[0], sizeof(Label_t *) * N / _K, cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();
  int expanded = 0;
  cuda_err_chk(cudaGetLastError());

  int DOP;
  std::vector<double> frontier_size;
  DOP = 100000;
  frontier_size.push_back(0);
  double *coefficients = (double *) malloc(AR_DEGREE * sizeof(double));
  int frontier_sum = 0;
  int max_bucket_size = 0;
  for (int i = 0; i < ntask; i++) {
    int total_label = 0;
    int MAX_LABELS_IN_FRONTIER = max_frontier(frontier_sum);

    DOP = MAX_LABELS_IN_FRONTIER / ITER;
    frontier_size.clear();
    int size_of_bucket = 1;
    int new_size_of_bucket = 0;
    cudaMemset(status, 0, sizeof(int) * N * SearchSpace);
    cudaMemset(bucket, 0, sizeof(Label_t) * N * SearchSpace);
    Label_t source = init_state_l[i];

    cudaMemcpy(bucket, &source, sizeof(Label_t), cudaMemcpyHostToDevice);
    cudaMemset(dominance_list, 0, sizeof(Label_t) * N * SearchSpace);
    cuda_err_chk(cudaGetLastError());

    cudaMemset(dynamic_list, 0, sizeof(Label_t) * N * SearchSpace * GlobalSpaceFactor);
    cudaMemset(dominance_list_len, 0, sizeof(int) * N);
    cudaDeviceSynchronize();
    cuda_err_chk(cudaGetLastError());
    cudaMemset(dynamic_list_len, 0, sizeof(int) * N / _K);
    ck_init<Label_t> <<< 1, 1 >>> (dominance_list_len, dominance_list_d, source);
    cudaDeviceSynchronize();
    cuda_err_chk(cudaGetLastError());
    thrust::device_ptr<int> status_ptr = thrust::device_pointer_cast(status);
    thrust::device_ptr<int> scan_out_ptr = thrust::device_pointer_cast(scan_out);
    frontier_sum = 0;
    double est = 0;
    bool f = false;
    int round = 0;
    cudaDeviceSynchronize();
    TimeMeasurer otherT;
    TimeMeasurer expandT;
    TimeMeasurer pruneT;
    TimeMeasurer t_ex;
    TimeMeasurer t_measure;

    t.tic();
    t_measure.resume();
    while (size_of_bucket != 0) {
      round++;
      cudaDeviceSynchronize();
      cuda_err_chk(cudaGetLastError());
      dim3 block(32, BLOCK_SIZE / 32);
      int num_of_block, SS;
      if (size_of_bucket > DOP) {
        expandT.resume();
//        printf(">>>>>>>>>> %d %d %d %lf\n", round, parallelism, size_of_bucket, parallelism * 1.0 / size_of_bucket);
        num_of_block = DOP / v_per_block;
        int left_frontier = size_of_bucket - DOP;
        if (DOP % v_per_block) {
          num_of_block++;
        }
        SS = N * SearchSpace / num_of_block;
        ck_expand_bucket_with_filter<Label_t, Edge_t> <<< num_of_block, block >>>(SS,
                                                                                  &bucket[left_frontier],
                                                                                  buffer,
                                                                                  status,
                                                                                  g->gg,
                                                                                  mutex,
                                                                                  dominance_list_d,
                                                                                  dominance_list_len,
                                                                                  C_l[0],
                                                                                  DOP,
                                                                                  dynamic_list_d,
                                                                                  dynamic_list_len);
        cudaDeviceSynchronize();
        cuda_err_chk(cudaGetLastError());
        expandT.pause();

        pruneT.resume();
        ck_recheck_bucket<Label_t, Edge_t><<<num_of_block, 128>>>(SS,                                                                           bucket,
                                                                           buffer,
                                                                           status,
                                                                           dominance_list_d,
                                                                           dominance_list_len,
                                                                           left_frontier,
                                                                           dynamic_list_d,
                                                                           dynamic_list_len);
        cudaDeviceSynchronize();
        frontier_sum += DOP;
        pruneT.pause();

      } else {
        expandT.resume();
        num_of_block = size_of_bucket / v_per_block;
        if (size_of_bucket % v_per_block) {
          num_of_block++;
        }
//        REACH_HERE;
//        printf("NUM blk %d %d %d\n", num_of_block, size_of_bucket, v_per_block);
        SS = N * SearchSpace / num_of_block;
        ck_expand_bucket_with_filter<Label_t, Edge_t> <<< num_of_block, block >>>(SS,
                                                                                  bucket,
                                                                                  buffer,
                                                                                  status,
                                                                                  g->gg,
                                                                                  mutex,
                                                                                  dominance_list_d,
                                                                                  dominance_list_len,
                                                                                  C_l[0],
                                                                                  size_of_bucket,
                                                                                  dynamic_list_d,
                                                                                  dynamic_list_len);
        cudaDeviceSynchronize();
        frontier_sum += size_of_bucket;
        expandT.pause();
      }
      cuda_err_chk(cudaGetLastError());
      otherT.resume();
      ////////// bucket Time
      thrust::inclusive_scan(status_ptr, status_ptr + num_of_block, scan_out_ptr);
      cudaDeviceSynchronize();
      cudaMemcpy(&new_size_of_bucket, &scan_out[num_of_block - 1], sizeof(int), cudaMemcpyDeviceToHost);
      cudaDeviceSynchronize();
      max_bucket_size = max(max_bucket_size, new_size_of_bucket);
      ////////// bucket Time
      otherT.pause();
//      printf("%d %d\n", round, new_size_of_bucket);

      // change K size
      double log_bucket = log2(new_size_of_bucket);
      if (log_bucket < frontier_size.back()
          && frontier_size.back() < frontier_size[frontier_size.size() - 2]
          && frontier_size[frontier_size.size() - 2] < frontier_size[frontier_size.size() - 3]) {
        // ending part
        DOP = min(MAX_LABELS_IN_FRONTIER, DOP + DOP / 2);
        frontier_size.push_back(log_bucket);
      } else {
        frontier_size.push_back(log_bucket);
        if (frontier_size.size() > AR_DEGREE) {
          if (!AR::AutoRegression(&frontier_size[0], frontier_size.size(), AR_DEGREE, coefficients)) {
            printf("AR error!\n");
          } else {
            est = 0;
            for (int j = 0; j < AR_DEGREE; j++) {
              est += coefficients[j] * frontier_size[frontier_size.size() - j - 1];
            }
            est = exp2(est);
//            printf("EST is %lf %d\n", est, new_size_of_bucket);
            if (est > new_size_of_bucket) {
              double p = calculate_p<Label_t>(est, new_size_of_bucket, MAX_LABELS_IN_FRONTIER);
              DOP = (int) std::round(p);
//              parallelism = max((int) (std::round(p / NUM_SM / v_per_block)), 1) * NUM_SM * v_per_block;
//              printf("EST is %lf %d\n", est, new_size_of_bucket);
//              printf("<><><>< %d %lf\n", parallelism, p);
              f = true;
            }
          }
        }
      }

      if (new_size_of_bucket == 0) {
        break;
      }
      otherT.resume();
      ////////// bucket Time
      ck_fill_buckets<Label_t> << < num_of_block, BLOCK_SIZE >> > (SS, buffer, status, scan_out, bucket);
      cudaDeviceSynchronize();
      size_of_bucket = new_size_of_bucket;
      ////////// bucket Time
      otherT.pause();
    }
    t.toc();
    t.print_ms("Execution Time");
    cudaDeviceSynchronize();

    // Get result back
//    int *dm_size = (int *) malloc(sizeof(int) * N);
//    cudaMemcpy(dm_size, dynamic_list_len, sizeof(int) * N / _K, cudaMemcpyDeviceToHost);

//    int pos = SearchSpace * sink_l[0];
//    auto t1 = dominance_list[pos].d;
//    auto t2 = dominance_list[pos].delay;
//    for (int i = 0; i < dominance_list_len[sink_l[0]]; i++) {
//      if (dominance_list[pos+i].d < t1) {
//        t1 = dominance_list[pos+i].d;
//        t2 = dominance_list[pos+i].delay;
//      }
//    }
//    std::cout << std::fixed << "<><><><><><><>" << t2 << " " << t1 << std::endl;


//    t.print_ns("Execution Time");
//    std::cout << "Bucket management: " << double (otherT.total_time_ns()) / 1000000 << " ms" << std::endl;
//    std::cout << "Expand Time: " << double (t_ex.total_time_ns()) / 1000000 << " ms" << std::endl;
//    std::cout << "Prune Time: " << double (pruneT.total_time_ns() + (expandT.total_time_ns() - t_ex.total_time_ns())) / 1000000 << " ms" << std::endl;
//    std::cout << "K == " << parallelism << std::endl;
    cudaDeviceSynchronize();
  }
  return 0;
}

template<typename Label_t>
__global__ void ck_CSPP_path(int source, Label_t **dominance_list, int *list_len, int *path, int *mutex) {
  __shared__ int len;
  len = *list_len;
  __shared__
  Label_t final_label;
  final_label.d = LRG_FLOAT;

  // find the best result in dominance_list
  __syncthreads();
  if (len == 0) {
//    printf("%d\n", source);
    return;
  }
  int x = threadIdx.x;
  __syncthreads();
  while (x < len) {
//    printf("%d -- %.2lf %d\n", x, (*dominance_list)[x].d, (*dominance_list)[x].delay);
    if (final_label.d > (*dominance_list)[x].d) {
      bool leave = true;
      while (leave) {
        if (atomicCAS(mutex, 0, 1) == 0) {
          if (final_label.d > (*dominance_list)[x].d) {
            Label_t::copy(&final_label, &(*dominance_list)[x]);
            __threadfence();
          }
          atomicExch(mutex, 0);
          leave = false;
        }
      }
    }
    x += blockDim.x;
  }
  __syncthreads();
  // generate path
  if (threadIdx.x == 0) {
    Label_t *ptr_label = &final_label;
    int h = 0;
    while (ptr_label != NULL) {
//      printf("%d ", ptr_label->father_id);
      path[h++] = ptr_label->father_id;
      ptr_label = ptr_label->pre;
    }
//    printf("%d\n", source);
  }
  __syncthreads();
}

#endif
