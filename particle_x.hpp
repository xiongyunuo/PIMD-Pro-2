#ifndef PARTICLE_X_HPP
#define PARTICLE_X_HPP

#include <limits>

#define D 2
#define M 8
#define kB 1
#define hBar 1

typedef double XNum;

struct XParticle {
  XNum coor[D];
  XNum vel[D];
  XNum m;
  int n, p;
  XNum theta[D][M];
  XNum vtheta[D][M];
  XNum Q[D][M];
};

inline bool IsNan(XNum val) {
  return val != val;
}

inline bool IsInf(XNum val) {
  XNum max = std::numeric_limits<XNum>::max();
  return !(-max <= val && val <= max);
}

#endif