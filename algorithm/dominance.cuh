#ifndef _dominance_cuh_
#define _dominance_cuh_

#include <cuda_runtime.h>
#include "../include/label.hpp"
#include "../include/constraint.hpp"
#include "../include/edge.hpp"
#include "../include/graph.hpp"
#include "../include/util.hpp"

#define SearchSpace 128
#define GlobalSpaceFactor 2
#define _K 4

template<typename Label_t>
__device__ void second_level_pruning(
    Label_t *cur_v,
    Label_t *v_dominance_list,
    int &dominance_list_len,
    Label_t *v_dynamic_list,
    int &dynamic_list_len,
    Label_t *local_buffer,
    int *expand_buffer_size,
    int *mutex) {
  __shared__ int flag;
  // -1 is dominated
  // 1 dominate others
  // 0 non, append

  flag = 0;
  int vid = cur_v->vid;
  // each thread now work on 1 comparison
  int lid = threadIdx.y * blockDim.x + threadIdx.x;
  if (dominance_list_len <= SearchSpace) {
    while (lid < dominance_list_len && lid < SearchSpace && flag != -1) {
      if (Label_t::dominate(v_dominance_list[lid], *cur_v)) {
        flag = -1;
      }
      lid += blockDim.x * blockDim.y;
    }
  } else {
    while (lid < SearchSpace && flag != -1) {
      if (Label_t::dominate(v_dominance_list[lid], *cur_v)) {
        flag = -1;
      }
      lid += blockDim.x * blockDim.y;
    }
    lid = threadIdx.y * blockDim.x + threadIdx.x;
    while (lid < dynamic_list_len && flag != -1) {
      if (Label_t::dominate_vid(v_dynamic_list[lid], *cur_v)) {
        flag = -1;
      }
      lid += blockDim.x * blockDim.y;
    }
  }
  __syncthreads();
  if (flag == -1) {
    cur_v->vid = -1;
    return;
  }
  lid = threadIdx.y * blockDim.x + threadIdx.x;
  if (dominance_list_len <= SearchSpace) {
    while (lid < dominance_list_len && lid <= SearchSpace && flag != 1) {
      // analysis dominance situation
      if (Label_t::dominate(*cur_v, v_dominance_list[lid])) {
        bool leave = true;
        while (leave) {
          if (atomicCAS(mutex + vid, 0, 1) == 0) {
            if (Label_t::dominate(*cur_v, v_dominance_list[lid])) {
              Label_t::copy(&v_dominance_list[lid], cur_v);
//            cur_v->pre = &v_dominance_list[lid];
              flag = 1;
            }
            __threadfence_block();
            atomicExch(mutex + vid, 0);
            leave = false;
          }
        }
      }
      lid += blockDim.x * blockDim.y;
    }
  } else {
    while (lid < SearchSpace && flag != 1) {
      // analysis dominance situation
      if (Label_t::dominate(*cur_v, v_dominance_list[lid])) {
        bool leave = true;
        while (leave) {
          if (atomicCAS(mutex + vid, 0, 1) == 0) {
            if (Label_t::dominate(*cur_v, v_dominance_list[lid])) {
              Label_t::copy(&v_dominance_list[lid], cur_v);
              flag = 1;
            }
            __threadfence_block();
            atomicExch(mutex + vid, 0);
            leave = false;
          }
        }
      }
      lid += blockDim.x * blockDim.y;
    }
    __syncthreads();
    lid = threadIdx.y * blockDim.x + threadIdx.x;
    while (lid < dynamic_list_len && flag != 1) {
      // analysis dominance situation
      if (Label_t::dominate_vid(*cur_v, v_dynamic_list[lid])) {
        bool leave = true;
        while (leave) {
          if (atomicCAS(mutex + (vid / _K), 0, 1) == 0) {
            if (Label_t::dominate_vid(*cur_v, v_dynamic_list[lid])) {
              Label_t::copy(&v_dynamic_list[lid], cur_v);
              flag = 1;
            }
            __threadfence_block();
            atomicExch(mutex + (vid / _K), 0);
            leave = false;
          }
        }
      }
      lid += blockDim.x * blockDim.y;
    }
  }
  __syncthreads();
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    if (flag == 0) {
      bool leave = true;
      while (leave) {
        if (dominance_list_len < SearchSpace) {
          if (atomicCAS(mutex + vid, 0, 1) == 0) {
            Label_t::copy(&v_dominance_list[dominance_list_len], cur_v);
            (dominance_list_len)++;
            __threadfence_block();
            atomicExch(mutex + vid, 0);
            leave = false;
          }
        } else {
          if (atomicCAS(mutex + (vid / _K), 0, 1) == 0) {
            Label_t::copy(&v_dynamic_list[dynamic_list_len], cur_v);
            (dynamic_list_len)++;
            (dominance_list_len)++;
            __threadfence_block();
            atomicExch(mutex + (vid / _K), 0);
            leave = false;
          }
        }
      }
    }
    // fill into bucket
    Label_t::copy(&local_buffer[*expand_buffer_size], cur_v);
    (*expand_buffer_size)++;
    cur_v->vid = -1;
  }
  __syncthreads();
}

template<typename Label_t>
__device__
void second_level_pruning_only(
    Label_t *cur_v,
    Label_t *v_dominance_list,
    int &dominance_list_len,
    Label_t *v_dynamic_list,
    int &dynamic_list_len,
    Label_t *local_buffer,
    int *expand_buffer_size) {
  __shared__ int flag;
  // -1 is dominated
  // 1 dominate others
  // 0 non, append
  flag = 0;
  int lid = threadIdx.y * blockDim.x + threadIdx.x;
  while (lid < dominance_list_len && lid < SearchSpace && flag != -1) {
    if (Label_t::dominate(v_dominance_list[lid], *cur_v) && !Label_t::equal(&v_dominance_list[lid], cur_v)) {
      flag = -1;
    }
    lid += blockDim.x * blockDim.y;
  }
  __syncthreads();
  if (flag == -1) {
    cur_v->vid = -1;
    return;
  }
  if (threadIdx.x == 0) {
    // fill into bucket
    Label_t::copy(&local_buffer[*expand_buffer_size], cur_v);
    (*expand_buffer_size)++;
    cur_v->vid = -1;
  }
  __syncthreads();
}

#endif
