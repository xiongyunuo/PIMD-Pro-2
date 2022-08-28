#ifndef RANDOM_X_HPP
#define RANDOM_X_HPP

typedef unsigned int XRandUInt;
typedef float XRandF;

void XSetRandSeed(XRandUInt seed);
XRandUInt XGetRandSeed();
XRandF XRandFloat();
XRandF XRandGauss();
XRandUInt XRandInteger();

#endif