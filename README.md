# Accelerating Exact Constrained Shortest Paths on GPUs

[![GitHub license](https://img.shields.io/badge/license-apache2-yellowgreen)](./LICENSE)

Vine is a a novel and practical framework for accelerating exact CSPs on the GPU.


Organization
--------

This repository contains code for our VLDB submission "Accelerating Exact Constrained Shortest Paths on GPUs".

- root folder contains main.cu and CMAkeLists.txt

- **include** folder contains helper functions, including graph storage, IO, CSP APIs, Timer, AR model, and others.

- **algorithm** folder contains CUDA and CPU-based CSP processing procedures.

- **inputs** sample input folders, the executable will read input/config for graph type, graph file, and query file.

- **script** folder contains python scripts that helps with reproduce the experiements

Compilation
--------

Requirements:
* nvcc > 7.0
* CMake >= 3.12


To build:
```sh
$ mkdir build
$ cd build && cmake ..
$ make
```

To run:

```
$ ./vine
```

The binary ```Vine``` will read the config file at ```./inputs/config``` and get the file path according to the config file. To run multiple query instances, please use the python scripts in the script folder.

Input Formats
-----------
We support the road network format used by [COLA](https://sourceforge.net/projects/cola2016/), and graphs from  [SNAP dataset collection](https://snap.stanford.edu/data/index.html). Please check main.cu and include/io.hpp for graph input.

Please check input/README.md for details about the input config file


Programming with Vine
---------------------
To use the APIs provided in ```Vine```, please refer to the ```label``` class defined in ```include/label.hpp```. Generally, a CSP solver can be programmed by providing t he following three functions, **expand**, **dominate**. The feasibility checking of the label is refactored to another class in ```include/contraint.hpp```. Please check the **pass_constraints_check** function as well.

```c++
class c_label {
 public:
  int vid; float d; float delay;
  __device__ c_label(int vid, float d, float delay) : vid(vid), d(d), delay(delay) {
    return;
  }

  template<typename Edge_t>
  __device__ void expand(int next_v, Edge_t e, c_label *son) {
    son->vid = next_v;
    son->d = this->d + e.len;
    son->delay = this->delay + e.delay;
  }

  __device__ static bool dominate(c_label a, c_label b) {
    return a.d <= b.d && a.delay <= b.delay;
  }
};
```



How to cite **Vine** 
-----------

If you use **Vine**, please cite our work ([full version](https://www.comp.nus.edu.sg/~hebs/pub/vine20.pdf)).
```
@article{shengliang2021vine,
  author = {Lu, Shengliang and He, Bingsheng and Li, Yuchen and Fu, Hao},
  title = {Accelerating Exact Constrained Shortest Paths on GPUs},
  journal = {PVLDB (Proceedings of the VLDB Endowment)},
  year = {2021}
}
```

## Related work
* Graph systems on GPU: [Medusa](https://github.com/Xtra-Computing/Medusa) | [Gunrock](https://github.com/gunrock/gunrock/) 
* Effective Indexing for Approximate Constrained Shortest Path Queries on Large Road Networks: [COLA technical report](https://sites.google.com/site/colatr2016/)
