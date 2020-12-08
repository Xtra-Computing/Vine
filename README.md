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

Compiler:
* nvcc > 7.0

Build system:
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

Input Formats
-----------
We support the road network format used by [COLA](https://sourceforge.net/projects/cola2016/), and graphs from  [SNAP dataset collection](https://snap.stanford.edu/data/index.html). Please check main.cu and include/io.hpp for graph input.

Please check input/README.md for details about the input config file

## How to cite **Vine** 
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
