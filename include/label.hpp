#ifndef _label_hpp_
#define _label_hpp_
#include <cuda_runtime.h>
#include "edge.hpp"
#include "io.hpp"
#ifndef CSP
#define CSP 1
#endif

__device__ float *A_dd, *B_dd, *Theta_dd;

__global__ void assign(float *A, float *B, float *Theta) {
  A_dd = A;
  B_dd = B;
  Theta_dd = Theta;
}

class label {
 public:
  int vid;
 public:
  __device__ label() {
    return;
  }

  __device__ virtual int get_vid() {
    return vid;
  }

  __device__ label
  expand(edge e) {
    return label();
  }

  __device__ static bool pass_constraints_check() {
    return false;
  }
  __device__ static bool dominate(label a, label b) {
    return false;
  }

  ~label() {
  }
};

#if CSP == 2
class c_label {
public	:
    int vid;
    float d;
    float delay;
    int hop;
    int father_id;
    c_label *pre;

public	:
    __device__ c_label() {
        return;
    }

    c_label(int vid) : vid(vid) {
        d = 0;
        delay = 0;
        hop = 0;
        father_id = -1;
        pre = NULL;
    }

    __device__ c_label(int vid, float d, float delay, int hop) :
    vid(vid), d(d), delay(delay), hop(hop) {
        pre = NULL;

        return;
    }

    __device__ static void copy(c_label *a, c_label *b) {
        a->vid = b->vid;
        a->d = b->d;
        a->delay = b->delay;
        a->hop = b->hop;
        a->father_id = b->father_id;
        a->pre = b->pre;
    }

    __device__ int get_vid() {
        return vid;
    }

    template<typename Edge_t>
    __device__ void expand(int next_v, Edge_t e, c_label *son) {
        son->vid = next_v;
        son->d = this->d + e.len;
        son->delay = this->delay + e.delay;
        son->hop = this->hop + 1;
        son->father_id = this->vid;
        son->pre = this->pre;
    }

    __device__ static bool dominate(c_label a, c_label b) {
        if (a.d <= b.d && a.delay <= b.delay && a.hop <= b.hop) {
            return true;
        }
        return false;
    }
};
#elif CSP == 3
class c_label {
public	:
    int vid;
    float d;
    float delay;
    int hop;
    float maxcost;
    int father_id;
    c_label *pre;

public	:
    __device__ c_label() {
        return;
    }

    c_label(int vid) : vid(vid) {
        d = 0;
        delay = 0;
        hop = 0;
        maxcost = 0;
        father_id = -1;
        pre = NULL;
    }

    __device__ c_label(int vid, float d, float delay, int hop, float maxcost) :
    vid(vid), d(d), delay(delay), hop(hop), maxcost(maxcost) {
        pre = NULL;

        return;
    }

    __device__ static void copy(c_label *a, c_label *b) {
        a->vid = b->vid;
        a->d = b->d;
        a->delay = b->delay;
        a->hop = b->hop;
        a->maxcost = b->maxcost;
        a->father_id = b->father_id;
        a->pre = b->pre;
    }

    __device__ int get_vid() {
        return vid;
    }

    template<typename Edge_t>
    __device__ void expand(int next_v, Edge_t e, c_label *son) {
        son->vid = next_v;
        son->d = this->d + e.len;
        son->delay = this->delay + e.delay;
        son->hop = this->hop + 1;
        son->maxcost = max(this->maxcost, e.len);
        son->father_id = this->vid;
        son->pre = this->pre;
    }

    __device__ static bool dominate(c_label a, c_label b) {
        if (a.d <= b.d && a.delay <= b.delay && a.hop <= b.hop && a.maxcost <= b.maxcost) {
            return true;
        }
        return false;
    }
};
#elif CSP == 1
class c_label {
 public    :
  int vid;
  float d;
  float delay;
//  int father_id;
//  c_label *pre;

 public    :
  __device__ c_label() {
    return;
  }

  c_label(int vid) : vid(vid) {
    d = 0;
    delay = 0;
//    father_id = -1;
//    pre = NULL;
  }

  __device__ c_label(int vid, float d, float delay) :
      vid(vid), d(d), delay(delay) {
//    pre = NULL;

    return;
  }
  __device__ static bool equal(c_label *a, c_label *b) {
    return a->vid == b->vid && a->d == b->d && a->delay == b->delay;
  }

  __device__ static void copy(c_label *a, c_label *b) {
    a->vid = b->vid;
    a->d = b->d;
    a->delay = b->delay;
//    a->father_id = b->father_id;
//    a->pre = b->pre;
  }

  __device__ int get_vid() {
    return vid;
  }

  template<typename Edge_t>
  __device__ void expand(int next_v, Edge_t e, c_label *son) {
    son->vid = next_v;
    son->d = this->d + e.len;
    son->delay = this->delay + e.delay;
//    son->father_id = this->vid;
//    son->pre = this->pre;
//    printf("%d->%d %lf %lf %lf %lf %lf %lf\n", this->vid, next_v, e.delay, this->delay, son->delay, this->d, e.len, son->d);
//    printf("%d->%d %lf %lf\n", this->vid, next_v, son->d, son->delay);
  }

  template<typename Edge_t>
  void expand_cpu(int next_v, Edge_t e, c_label *son) {
    son->vid = next_v;
    son->d = this->d * e.len;
    son->delay = this->delay * e.delay;
//    son->father_id = this->vid;
//    son->pre = this->pre;
//    printf("%d->%d %lf %lf %lf\n", this->vid + 1, next_v + 1, e.delay, this->delay, son->delay);
  }

  __device__ static bool dominate(c_label a, c_label b) {
    return a.d <= b.d && a.delay <= b.delay;
  }
  __device__ static bool dominate_vid(c_label a, c_label b) {
    if (a.vid != b.vid) {
      return false;
    }
    return a.d <= b.d && a.delay <= b.delay;
  }

  static bool dominate_cpu(c_label a, c_label b) {
    return a.d <= b.d && a.delay <= b.delay;
  }
};
#elif CSP == 4 //2-OSP
class c_label {
public	:
    int vid;
    float d;
    float delay;
    int father_id;
    c_label *pre;

public	:
    __device__ c_label() {
        return;
    }

    c_label(int vid) : vid(vid) {
        d = 0;
        delay = 0;
        father_id = -1;
        pre = NULL;
    }

    __device__ c_label(int vid, float d, float delay) :
    vid(vid), d(d), delay(delay) {
        pre = NULL;

        return;
    }

    __device__ static void copy(c_label *a, c_label *b) {
        a->vid = b->vid;
        a->d = b->d;
        a->delay = b->delay;
        a->father_id = b->father_id;
        a->pre = b->pre;
    }

    __device__ int get_vid() {
        return vid;
    }

    template<typename Edge_t>
    __device__ void expand(int next_v, Edge_t e, c_label *son) {
        son->vid = next_v;
        son->d = this->d + e.len;
        son->delay = this->delay + e.delay;
        son->father_id = this->vid;
        son->pre = this->pre;
    }

    __device__ static bool dominate(c_label a, c_label b) {
        if (a.d <= b.d && a.delay <= b.delay) {
            return true;
        }
        return false;
    }
};
#elif CSP == 5 // social
class c_label {
 public    :
  int vid;
  float d;
  float delay;
//  int father_id;
//  c_label *pre;

 public    :
  __device__ c_label() {
    return;
  }

  c_label(int vid) : vid(vid) {
    d = 1;
    delay = 1;
//    father_id = -1;
//    pre = NULL;
  }

  __device__ c_label(int vid, float d, float delay) :
      vid(vid), d(d), delay(delay) {
//    pre = NULL;

    return;
  }
  __device__ static bool equal(c_label *a, c_label *b) {
    return a->vid == b->vid && a->d == b->d && a->delay == b->delay;
  }

  __device__ static void copy(c_label *a, c_label *b) {
    a->vid = b->vid;
    a->d = b->d;
    a->delay = b->delay;
//    a->father_id = b->father_id;
//    a->pre = b->pre;
  }

  __device__ int get_vid() {
    return vid;
  }

  template<typename Edge_t>
  __device__ void expand(int next_v, Edge_t e, c_label *son) {
    son->vid = next_v;
    son->d = this->d * e.len;
    son->delay = this->delay * e.delay;
//    son->father_id = this->vid;
//    son->pre = this->pre;
//    printf("%d->%d %lf %lf %lf %lf %lf %lf\n", this->vid, next_v, e.delay, this->delay, son->delay, this->d, e.len, son->d);
//    printf("%d->%d %lf %lf\n", this->vid, next_v, son->d, son->delay);
  }

  template<typename Edge_t>
  void expand_cpu(int next_v, Edge_t e, c_label *son) {
    son->vid = next_v;
    son->d = this->d * e.len;
    son->delay = this->delay * e.delay;
//    son->father_id = this->vid;
//    son->pre = this->pre;
//    printf("%d->%d %lf %lf %lf\n", this->vid + 1, next_v + 1, e.delay, this->delay, son->delay);
  }

  __device__ static bool dominate(c_label a, c_label b) {
    return a.d >= b.d && a.delay >= b.delay;
  }
  __device__ static bool dominate_vid(c_label a, c_label b) {
    if (a.vid != b.vid) {
      return false;
    }
    return a.d >= b.d && a.delay >= b.delay;
  }

  static bool dominate_cpu(c_label a, c_label b) {
    return a.d >= b.d && a.delay >= b.delay;
  }
};
#else // network qos
class c_label {
 public    :
  int vid;
  float d;
  float delay;
//  int father_id;
//  c_label *pre;

 public    :
  __device__ c_label() {
    return;
  }

  c_label(int vid) : vid(vid) {
    d = 0;
    delay = 21474836;
//    father_id = -1;
//    pre = NULL;
  }

  __device__ c_label(int vid, float d, float delay) :
      vid(vid), d(d), delay(delay) {
//    pre = NULL;

    return;
  }
  __device__ static bool equal(c_label *a, c_label *b) {
    return a->vid == b->vid && a->d == b->d && a->delay == b->delay;
  }

  __device__ static void copy(c_label *a, c_label *b) {
    a->vid = b->vid;
    a->d = b->d;
    a->delay = b->delay;
//    a->father_id = b->father_id;
//    a->pre = b->pre;
  }
  __device__ int get_vid() {
    return vid;
  }

  template<typename Edge_t>
  __device__ void expand(int next_v, Edge_t e, c_label *son) {
    son->vid = next_v;
    son->d = this->d + e.len;
    son->delay = min(this->delay, e.delay);
//    son->father_id = this->vid;
//    son->pre = this->pre;
//    printf("%d->%d %lf %lf %lf %lf %lf %lf\n", this->vid, next_v, e.delay, this->delay, son->delay, this->d, e.len, son->d);
//    printf("%d->%d %lf %lf\n", this->vid, next_v, son->d, son->delay);
  }

  template<typename Edge_t>
  void expand_cpu(int next_v, Edge_t e, c_label *son) {
    son->vid = next_v;
    son->d = this->d + e.len;
    son->delay = min(this->delay, e.delay);
//    son->father_id = this->vid;
//    son->pre = this->pre;
//    printf("%d->%d %lf %lf %lf\n", this->vid + 1, next_v + 1, e.delay, this->delay, son->delay);
  }

  __device__ static bool dominate(c_label a, c_label b) {
    if (a.d <= b.d && a.delay >= b.delay) {
        return true;
    }
    return false;
  }
  __device__ static bool dominate_vid(c_label a, c_label b) {
    if (a.vid != b.vid) {
      return false;
    }
    if (a.d <= b.d && a.delay >= b.delay) {
      return true;
    }
    return false;
  }

  static bool dominate_cpu(c_label a, c_label b) {
    return a.d <= b.d && a.delay >= b.delay;
  }
};
#endif

#endif