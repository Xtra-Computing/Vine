file(GLOB_RECURSE CURRENT_HEADERS *.hpp)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_60 -std=c++11 -Xptxas -dlcm=cg)
cuda_add_library(include_lib ../empty.cpp AR.cpp ${CURRENT_HEADERS})