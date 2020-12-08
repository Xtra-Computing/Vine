#ifndef _bucket_cuh_
#define _bucket_cuh_

#include <cuda_runtime.h>
#include "dominance.cuh"
#include "../include/gpu_graph.hpp"

#define LRG_FLOAT 1000000
#define BLOCK_SIZE 32
#define TABLE_SIZE 96
#define VSIZE 2

#define HASH

template<typename Label_t, typename Edge_t>
__global__ void ck_expand_bucket_with_filter(
    int local_size,
    Label_t *bucket,
    Label_t *buffer,
    int *st,
    gpuGraph<Edge_t> *gg,
    int *mutex,
    Label_t **dominance_list,
    int *dominance_list_len,
    c_constraint<Label_t> C,
    int size_of_bucket,
    Label_t **dynamic_list,
    int *dynamic_list_len) {
  __shared__ Label_t hash_table[TABLE_SIZE];
  __shared__ int maxR;
#ifdef HASH
  __shared__ int mtx;
  __shared__ int hs_cout;
  mtx = 0;
  hs_cout = 0;
#endif
  int tid = threadIdx.y * blockDim.x + threadIdx.x;
  while (tid < TABLE_SIZE) {
    hash_table[tid].vid = -1;
    tid += BLOCK_SIZE;
  }
  st[blockIdx.x] = 0;
  maxR = 0;
  int l, r;
  Label_t *local_buffer = buffer + blockIdx.x * local_size;
  Label_t *cur_v;
  tid = (blockIdx.x * blockDim.y + threadIdx.y) * VSIZE;
  for (int T = 0; T < VSIZE; T++, tid++) {
    if (tid < size_of_bucket) {
      cur_v = bucket + tid;
      l = gg->_count_row_d[cur_v->vid];
      r = gg->_count_row_d[cur_v->vid + 1];
    } else {
      cur_v = nullptr;
      l = r = -1; // give up processing the last empty frontier
    }
    // get the max rounds
    __syncthreads();
    if (threadIdx.x == 0) {
      for (int i = 0; i < BLOCK_SIZE / 32; i++) {
        if (r - l > maxR) {
          maxR = r - l;
        }
      }
    }
    __syncthreads();
//    if (threadIdx.x % 32 == 0) {
//      printf("id %d %d %d %d %d %d\n", blockIdx.x, threadIdx.y * 32 + threadIdx.x, blockIdx.x * blockDim.y + threadIdx.y, tid, r, maxR);
//    }
    l = threadIdx.x + l;
    for (; maxR > 0; maxR -= blockDim.x) {
      if (l < r) {
        Edge_t *next_e = gg->el_d + gg->_eid_col_d[l];
        Label_t new_label;
        cur_v->expand(gg->_vid_col_d[l], *next_e, &new_label);
        if (!C.pass_constraints_check(new_label)) {
          break;
        }
#ifdef HASH
        int pos = new_label.vid % TABLE_SIZE;

//          printf("id %d %d %d ***\n", tid, new_label.vid, pos);
        bool inserted = false;
        while (!inserted) {
          if (hash_table[pos].vid == -1) {
            bool leave = true;
            while (leave) {
              if (atomicCAS(&mtx, 0, 1) == 0) {
                if (hash_table[pos].vid == -1) {
                  Label_t::copy(&hash_table[pos], &new_label);
                  hs_cout++;
                  inserted = true;
                }
                __threadfence();
                atomicExch(&mtx, 0);
                leave = false;
              }
            }
          }
          if (inserted || Label_t::dominate_vid(hash_table[pos], new_label)) {
            break;
          }
          if (Label_t::dominate_vid(new_label, hash_table[pos])) {
            bool leave = true;
            while (leave) {
              if (atomicCAS(&mtx, 0, 1) == 0) {
                if (Label_t::dominate_vid(new_label, hash_table[pos])) {
                  Label_t::copy(&hash_table[pos], &new_label);
                  hs_cout++;
                  inserted = true;
                }
                __threadfence();
                atomicExch(&mtx, 0);
                leave = false;
              }
            }
          }
          pos = (pos + 1) % TABLE_SIZE;
        }
#else
        hash_table[threadIdx.y * blockDim.x + threadIdx.x] = new_label;
        for (int j = 0; j < blockDim.x * blockDim.y; j++) {
            if (hash_table[j].vid > 0) {
                int vid = hash_table[j].vid;
                second_level_pruning<Label_t>(&hash_table[j],
                                              dominance_list[vid],
                                              dominance_list_len[vid],
                                              dynamic_list[vid / _K],
                                              dynamic_list_len[vid / _K],
                                              local_buffer,
                                              &st[blockIdx.x],
                                              mutex);
            }
        }
#endif
      }
      __syncthreads();
#ifdef HASH
      if (hs_cout + BLOCK_SIZE > TABLE_SIZE) {
        hs_cout = 0;
#ifdef DUPLICATE_EDGE
        // works only for duplicated edge
       for (int j = 0; j < BLOCK_SIZE; j++) {
         if (new_label[j].vid < 0) {
           break;
         }
         if (threadIdx.x != j) {
           if (new_label[threadIdx.x].vid == new_label[j].vid && Label_t::dominate(new_label[threadIdx.x], new_label[j])) {
             new_label[j].delay = LRG_FLOAT;
           }
         }
         __syncthreads();
       }
#endif
        for (int j = 0; j < TABLE_SIZE; j++) {
          int vid = hash_table[j].vid;
          if (vid >= 0) {
            second_level_pruning<Label_t>(&hash_table[j],
                                          dominance_list[vid],
                                          dominance_list_len[vid],
                                          dynamic_list[vid / _K],
                                          dynamic_list_len[vid / _K],
                                          local_buffer,
                                          &st[blockIdx.x],
                                          mutex);

          }
        }
      }
#else
      ;;
#endif
      l += blockDim.x;
      __syncthreads();
    }
  }

#ifdef HASH

  for (int j = 0; j < TABLE_SIZE; j++) {
    int vid = hash_table[j].vid;
    if (vid >= 0) {
      second_level_pruning<Label_t>(&hash_table[j],
                                    dominance_list[vid],
                                    dominance_list_len[vid],
                                    dynamic_list[vid / _K],
                                    dynamic_list_len[vid / _K],
                                    local_buffer,
                                    &st[blockIdx.x],
                                    mutex);

    }
  }
#endif
}

template<typename Label_t>
__global__ void ck_fill_buckets(int local_size, Label_t *buffer, int *st, int *histogram, Label_t *bucket) {
  int base;
  int x, index;
  if (blockIdx.x == 0) {
    base = 0;
  } else {
    base = histogram[blockIdx.x - 1];
  }
  __syncthreads();
  index = blockIdx.x * local_size;
  x = threadIdx.x;
  __syncthreads();
  while (x < st[blockIdx.x]) {
    bucket[base + x] = buffer[index + x];
//    printf("%d %.3lf %.3lf\n", bucket[base + x].vid, bucket[base + x].d, bucket[base + x].delay);
    x += blockDim.x;
  }
}

template<typename Label_t, typename Edge_t>
__global__ void ck_recheck_bucket(
    int local_size,
    Label_t *bucket,
    Label_t *buffer,
    int *st,
    Label_t **dominance_list,
    int *dominance_list_len,
    int size_of_bucket,
    Label_t **dynamic_list,
    int *dynamic_list_len) {
  int tid = blockIdx.x;
  Label_t *local_buffer = buffer + blockIdx.x * local_size;
  while (tid < size_of_bucket) {
    int vid = bucket[tid].vid;
//        if (threadIdx.x == 0) {
//            printf("id %d %d %d\n", blockIdx.x, tid, vid);
//        }

    second_level_pruning_only<Label_t>(&bucket[tid],
                                       dominance_list[vid],
                                       dominance_list_len[vid],
                                       dynamic_list[vid / _K],
                                       dynamic_list_len[vid / _K],
                                       local_buffer,
                                       &st[blockIdx.x]);
    tid += gridDim.x;
  }
}

template<typename Label_t, typename Edge_t>
__global__ void ck_expand_bucket_with_nothing(
    int local_size,
    Label_t *bucket,
    Label_t *buffer,
    int *st,
    gpuGraph<Edge_t> *gg,
    int *mutex,
    Label_t **dominance_list,
    int *dominance_list_len,
    c_constraint<Label_t> C,
    int size_of_bucket,
    Label_t **dynamic_list,
    int *dynamic_list_len) {
  __shared__ Label_t hash_table[TABLE_SIZE];
  __shared__ int maxR;
#ifdef HASH
  __shared__ int mtx;
  __shared__ int hs_cout;
  mtx = 0;
  hs_cout = 0;
#endif
  int tid = threadIdx.y * blockDim.x + threadIdx.x;
  while (tid < TABLE_SIZE) {
    hash_table[tid].vid = -1;
    tid += BLOCK_SIZE;
  }
//  st[blockIdx.x] = 0;
  maxR = 0;
  int l, r;
//  Label_t *local_buffer = buffer + blockIdx.x * local_size;
  Label_t *cur_v;
  tid = (blockIdx.x * blockDim.y + threadIdx.y) * VSIZE;
  for (int T = 0; T < VSIZE; T++, tid++) {
    if (tid < size_of_bucket) {
      cur_v = bucket + tid;
      l = gg->_count_row_d[cur_v->vid];
      r = gg->_count_row_d[cur_v->vid + 1];
    } else {
      cur_v = nullptr;
      l = r = -1; // give up processing the last empty frontier
    }
    // get the max rounds
    __syncthreads();
    if (threadIdx.x == 0) {
      for (int i = 0; i < BLOCK_SIZE / 32; i++) {
        if (r - l > maxR) {
          maxR = r - l;
        }
      }
    }
    __syncthreads();
//        if (threadIdx.x % 32 == 0) {
//            printf("id %d %d %d %d %d %d\n", blockIdx.x, threadIdx.y * 32 + threadIdx.x,
//                   blockIdx.x * blockDim.y + threadIdx.y, tid, r, maxR);
//        }
    l = threadIdx.x + l;
    for (; maxR > 0; maxR -= blockDim.x) {
      if (l < r) {
        Edge_t *next_e = gg->el_d + gg->_eid_col_d[l];
        Label_t new_label;
        cur_v->expand(gg->_vid_col_d[l], *next_e, &new_label);
//        if (!C.pass_constraints_check(new_label)) {
//          break;
//        }
#ifdef HASH
        int pos = new_label.vid % TABLE_SIZE;

//          printf("id %d %d %d ***\n", tid, new_label.vid, pos);
        bool inserted = false;
        while (!inserted) {
          if (hash_table[pos].vid == -1) {
            bool leave = true;
            while (leave) {
              if (atomicCAS(&mtx, 0, 1) == 0) {
                if (hash_table[pos].vid == -1) {
//                  Label_t::copy(&hash_table[pos], &new_label);
                  hs_cout++;
                  inserted = true;
                }
                __threadfence();
                atomicExch(&mtx, 0);
                leave = false;
              }
            }
          }
          if (inserted || Label_t::dominate_vid(hash_table[pos], new_label)) {
            break;
          }
          if (Label_t::dominate_vid(new_label, hash_table[pos])) {
            bool leave = true;
            while (leave) {
              if (atomicCAS(&mtx, 0, 1) == 0) {
                if (Label_t::dominate_vid(new_label, hash_table[pos])) {
//                  Label_t::copy(&hash_table[pos], &new_label);
                  hs_cout++;
                  inserted = true;
                }
                __threadfence();
                atomicExch(&mtx, 0);
                leave = false;
              }
            }
          }
          pos = (pos + 1) % TABLE_SIZE;
        }
#else
        hash_table[threadIdx.y * blockDim.x + threadIdx.x] = new_label;
        for (int j = 0; j < blockDim.x * blockDim.y; j++) {
            if (hash_table[j].vid > 0) {
                int vid = hash_table[j].vid;
                //            second_level_pruning<Label_t>(&hash_table[j],
                //                                          dominance_list[vid],
                //                                          dominance_list_len[vid],
                //                                          dynamic_list[vid / _K],
                //                                          dynamic_list_len[vid / _K],
                //                                          local_buffer,
                //                                          &st[blockIdx.x],
                //                                          mutex);
            }
        }
#endif
      }
      __syncthreads();
#ifdef HASH
      if (hs_cout + BLOCK_SIZE > TABLE_SIZE) {
        hs_cout = 0;
#ifdef DUPLICATE_EDGE
        // works only for duplicated edge
       for (int j = 0; j < BLOCK_SIZE; j++) {
         if (new_label[j].vid < 0) {
           break;
         }
         if (threadIdx.x != j) {
           if (new_label[threadIdx.x].vid == new_label[j].vid && Label_t::dominate(new_label[threadIdx.x], new_label[j])) {
             new_label[j].delay = LRG_FLOAT;
           }
         }
         __syncthreads();
       }
#endif
        for (int j = 0; j < TABLE_SIZE; j++) {
          int vid = hash_table[j].vid;
          if (vid >= 0) {
//            second_level_pruning<Label_t>(&hash_table[j],
//                                          dominance_list[vid],
//                                          dominance_list_len[vid],
//                                          dynamic_list[vid / _K],
//                                          dynamic_list_len[vid / _K],
//                                          local_buffer,
//                                          &st[blockIdx.x],
//                                          mutex);

          }
        }
      }
#else
      ;;
#endif
      l += blockDim.x;
      __syncthreads();
    }
  }

#ifdef HASH

  for (int j = 0; j < TABLE_SIZE; j++) {
    int vid = hash_table[j].vid;
    if (vid >= 0) {
//      second_level_pruning<Label_t>(&hash_table[j],
//                                    dominance_list[vid],
//                                    dominance_list_len[vid],
//                                    dynamic_list[vid / _K],
//                                    dynamic_list_len[vid / _K],
//                                    local_buffer,
//                                    &st[blockIdx.x],
//                                    mutex);

    }
  }
#endif
}

#endif
