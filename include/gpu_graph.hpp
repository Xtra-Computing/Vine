#ifndef _gpu_graph_hpp_
#define _gpu_graph_hpp_
#include <cuda_runtime.h>

template<typename Edge_t>
class gpuGraph {
 public:
  int N;
  bool direct;

  int *_count_row_d; // CSR pointer
  int *_eid_col_d; // CSR edge id
  int *_vid_col_d; //CSR tail v id
  int *_sink_d; // eid to sink
  int *_csr_d; // eid to csr pos
  float *_len_data_d; // CSR edge len
  float *_delay_data_d; // CSR edge delay

  Edge_t *el_d;
 public:
  gpuGraph() {
    // cudaMalloc(gg_d);
    // cudaMemcpy(gg_d, g->gg);
  }
  ~gpuGraph() {
    //	cudaFree(csr_count_row_d);
    //	cudaFree(csr_eid_col_d);
    //	cudaFree(csr_vid_col_d);
    //	cudaFree(eid_cap_d);
    //	cudaFree(eid_sink_d);
    //	cudaFree(eid_csr_d);
    //	cudaFree(csr_len_data_d);
    //	cudaFree(csr_delay_data_d);
  }
};

#endif /* GPUGRAPH_H_ */
