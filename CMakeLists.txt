cmake_minimum_required(VERSION 3.12)
project(Vine)
set(USE_CUDA ON CACHE BOOL "Compile with CUDA")
find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_CXX_STANDARD 11)
set(CUDA_SEPARABLE_COMPILATION ON)
find_package(CUDA REQUIRED QUIET)
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_60 -std=c++11 -Xptxas -dlcm=cg)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_60 -G -g -std=c++11 -Xptxas -dlcm=cg)

add_definitions(-DCSP=1)
add_subdirectory(algorithm)
add_subdirectory(include)

cuda_add_executable(vine
        main.cu)
target_link_libraries(vine
        algorithm_lib
        include_lib)