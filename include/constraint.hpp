#ifndef _csp_hpp_
#define _csp_hpp_
#include <cuda_runtime.h>
#include "label.hpp"

template<typename Label_t>
class constraint {
 public:
  constraint() {
    return;
  }
  virtual void set_constraints() {
    return;
  }

  __device__ virtual bool pass_constraints_check(Label_t l) {
    return false;
  }
};

#if CSP == 2
template<typename Label_t>
class c_constraint {
public	:
    float delay_L;
    int hop_L;
public	:
    c_constraint() {
        return;
    }
    c_constraint(float d, int h = 100) :
    delay_L(d), hop_L(h) {
        return;
    }
    void set_constraints(float d, int h) {
        delay_L = d;
        hop_L = h;

        return;
    }

    __device__ bool pass_constraints_check(const Label_t &l) {
        if (l.delay > delay_L || l.hop > hop_L) {
            return false;
        }

        return true;
    }
};
#elif CSP == 3
template<typename Label_t>
class c_constraint {
public	:
    float delay_L;
    int hop_L;
    float maxcost_L;
public	:
    c_constraint() {
        return;
    }
    c_constraint(float d, int h = 100, float m = 10000) :
    delay_L(d), hop_L(h), maxcost_L(m) {
        return;
    }
    void set_constraints(float d, int h, int m) {
        delay_L = d;
        hop_L = h;
        maxcost_L = m;

        return;
    }

    __device__ bool pass_constraints_check(const Label_t &l) {
        if (l.delay > delay_L || l.hop > hop_L || l.maxcost > maxcost_L) {
            return false;
        }

        return true;
    }
};
#elif CSP == 1
template<typename Label_t>
class c_constraint {
 public    :
  float delay_L;
 public    :
  c_constraint() {
    return;
  }
  c_constraint(float d) : delay_L(d) {
    return;
  }
  void set_constraints(int d) {
    delay_L = d;

    return;
  }

  __device__ bool pass_constraints_check(const Label_t &l) {
    return !(l.delay > delay_L);
  }
  bool pass_constraints_check_cpu(const Label_t &l) {
    return !(l.delay > delay_L);
  }
};
#elif CSP == 4 // 2-OSP
template<typename Label_t>
class c_constraint {
public	:
    float delay_L;
public	:
    c_constraint() {
        return;
    }
    c_constraint(float d) :
    delay_L(d) {
        return;
    }
    void set_constraints(int d) {
        delay_L = d;

        return;
    }

    __device__ bool pass_constraints_check(const Label_t &l) {
        if (l.delay > delay_L) {
            return false;
        }

        return true;
    }
};

#elif CSP == 5 //Social
template<typename Label_t>
class c_constraint {
public	:
    float delay_L;
public	:
    c_constraint() {
        return;
    }
    c_constraint(float d) :
    delay_L(d) {
        return;
    }
    void set_constraints(int d) {
        delay_L = d;

        return;
    }

    __device__ bool pass_constraints_check(const Label_t &l) {
        if (l.delay > delay_L) {
            return false;
        }

        return true;
    }
};

#else // network qos
template<typename Label_t>
class c_constraint {
 public    :
  float delay_L;
 public    :
  c_constraint() {
    return;
  }
  c_constraint(float d) : delay_L(d) {
    return;
  }
  void set_constraints(int d) {
    delay_L = d;

    return;
  }

  __device__ bool pass_constraints_check(const Label_t &l) {
    return !(l.delay < delay_L);
  }
  bool pass_constraints_check_cpu(const Label_t &l) {
    return !(l.delay < delay_L);
 }
};
#endif

#endif
