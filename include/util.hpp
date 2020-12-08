#ifndef _util_hpp_
#define _util_hpp_

#include <iostream>

#define __DEBUG

#ifdef __DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
  }
}

#define cuda_err_chk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#endif

