#ifndef CSP_AR_HPP
#define CSP_AR_HPP
#define MAXENTROPY      0
#define LEASTSQUARES    1
#define TRUE  1
#define FALSE 0
#define MAXCOEFF 100

class AR {
 public:
  static int AutoRegression(double *, int, int, double *, int method = 0);
  static int ARMaxEntropy(double *, int, int, double **, double *, double *, double *, double *);
  static int ARLeastSquare(double *, int, int, double *);
  static int SolveLE(double **, double *, unsigned int);
};

#endif //CSP_AR_HPP
