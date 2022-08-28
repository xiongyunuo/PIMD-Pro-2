#ifndef STATISTIC_X_HPP
#define STATISTIC_X_HPP

#include "particle_x.hpp"
#include <complex>

class XStatistic {
public:
  XNum e;
  XNum T;
  std::complex<double> sign;
  std::complex<double> es;
  XNum *den;
  int num;
  XStatistic();
  void Init(int size);
  ~XStatistic() { delete[] den; }
};

#endif