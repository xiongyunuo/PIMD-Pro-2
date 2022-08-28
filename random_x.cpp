#include "random_x.hpp"
#include <cmath>

static XRandUInt XSeed = 1;
static XRandUInt XMod = 2147483647;
static XRandUInt XA = 16807;
static XRandUInt XQ = 127766;
static XRandUInt XR = 120485;

XRandUInt XRandInteger() {
  int res = (int)(XA * (XSeed % XQ)) - (int)(XR * (XSeed / XQ));
  if (res < 0) res += XMod;
  XSeed = res;
  return XSeed;
}

void XSetRandSeed(XRandUInt seed) {
  XSeed = seed;
}

XRandUInt XGetRandSeed() {
  return XSeed;
}

XRandF XRandFloat() {
  return (XRandF)(XRandInteger()) / XMod;
}

XRandF XRandGauss() {
  static XRandF first, second;
  static bool has = false;
  if (has) {
    has = false;
    return second;
  }
  has = true;
  XRandF r = std::sqrt(-2 * std::log(XRandFloat()));
  XRandF angle = 2 * M_PI * XRandFloat();
  first = r * cos(angle);
  second = r * sin(angle);
  return first;
}