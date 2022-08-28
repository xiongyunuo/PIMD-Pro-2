#include "statistic_x.hpp"

XStatistic::XStatistic() {
  e = 0;
  T = 0;
  sign = 0;
  es = 0;
}

void XStatistic::Init(int size) {
  int i;
  num = size;
  den = new XNum[size];
  for (i = 0; i < size; ++i)
    den[i] = 0;
}