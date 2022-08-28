#ifndef SIMULATION_X_HPP
#define SIMULATION_X_HPP

#include "particle_x.hpp"
#include "statistic_x.hpp"
#include <cmath>
#include <iostream>
#include <string>
#include <complex>

class XSimulation {
private:
  XParticle *particles;
  XParticle center;
  int N, P;
  XNum L;
  XNum beta, T;
  XNum vi;
  XNum beta2, vi2, g2, s2;
  XNum omgP;
  XNum omg0;
  XNum g, s;
  XNum **ENkCache;
  XNum *VBCache;
  std::complex<double> *VBCache2;
  XNum *ForceVBCache;
  XNum *ForceCache;
  XNum t, h;
  int count;
  XNum incre;
  XStatistic stat;
public:
  int skip, step;
  bool ok;
  inline XNum MinimumImage1(XNum a);
  inline XNum MinimumImage2(XNum a);
  inline int Index(int l, int j) { return (l-1)*P+j-1; }
  XNum Distance(XParticle *p1, XParticle *p2);
  void RelativeDistance(XParticle *p1, XParticle *p2, XNum *displace);
  XSimulation() {}
  void Initial(std::istream &in);
  void Dump(std::ostream &out);
  void VelocityRescale();
  XNum *Force();
  void UpdateMNHC_VV3();
  XNum ENk(int N2, int k);
  int NextIndex(int l, int j, int N2, int k);
  int PrevIndex(int l, int j, int N2, int k);
  XNum XExp(XNum k, XNum E, XNum EE);
  std::complex<double> XExp2(XNum k, std::complex<double> E, std::complex<double> EE);
  XNum XMinE(int N2);
  void FillVB();
  std::complex<double> XMinE2(int N2);
  void FillVB2();
  void dENk(int N2, int k, int l, int j, XNum *res);
  void FillForceVB();
  void FillENk();
  void TrapForce(int index, XNum *res);
  void PairForce(int index, int index2, XNum *res);
  void DumpForce(std::ostream &out);
  XNum VBEnergy();
  XNum TrapEnergy(int index);
  XNum PairEnergy(int index, int index2);
  XNum TotalEnergy();
  void DumpEnergy(std::ostream &out);
  void NHForce(XNum m, int f, XNum *v, XNum *Q, XNum *vtheta, XNum *res);
  void PeriodBoundary();
  XNum Temperature();
  void MakeStat();
  void DumpStat(std::ostream &out);
  void NormalizeStat();
  void SetOmg(XNum o) { omg0 = o; }
  XNum Partition();
  std::complex<double> Partition2();
  ~XSimulation();
};

XNum XSimulation::MinimumImage1(XNum a) {
  if (std::abs(a) > L / 2)
    return L - std::abs(a);
  return a;
}

XNum XSimulation::MinimumImage2(XNum a) {
  if (std::abs(a) > L / 2) {
    if (a < 0)
      return L - std::abs(a);
    else
      return -(L - std::abs(a));
  }
  return a;
}

#endif
