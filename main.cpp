#include "simulation_x.hpp"
#include <fstream>
#include <iostream>
#include <ctime>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Wrong number of arguments supplied" << std::endl;
    return 0;
  }
  std::ifstream in_file;
  in_file.open(std::string(argv[1]));
  std::ofstream out_file;
  out_file.open(std::string(argv[2]));
  XSimulation simul;
  /*simul.Initial(in_file);
  simul.Dump(out_file);
  simul.FillENk();
  simul.FillVB();
  simul.FillForceVB();
  simul.DumpForce(out_file);
  simul.DumpEnergy(out_file);
  simul.UpdateMNHC_VV3();
  simul.UpdateMNHC_VV3();
  //simul.Dump(out_file);
  for (int i = 0; i < 1000; ++i)
    simul.UpdateMNHC_VV3();
  simul.Dump(out_file);*/
  simul.Initial(in_file);
  int i;
  clock_t t;
  t = clock();
  //simul.Dump(out_file);
  for (i = 0; i < simul.skip; ++i) {
    simul.PeriodBoundary();
    simul.UpdateMNHC_VV3();
    if (!simul.ok) {
      simul.Dump(out_file);
      simul.DumpForce(out_file);
      simul.DumpEnergy(out_file);
      simul.NormalizeStat();
      simul.DumpStat(out_file);
      return 0;
    }
  }
  /*int n = simul.step / 5;
  for (i = 1; i <= simul.step; ++i) {
    simul.PeriodBoundary();
    simul.UpdateMNHC_VV3();
    if (!simul.ok) {
      simul.Dump(out_file);
      simul.DumpForce(out_file);
      simul.DumpEnergy(out_file);
      simul.NormalizeStat();
      simul.DumpStat(out_file);
      return 0;
    }
    if (i % 100000 == 0)
      out_file << "i=" << i << std::endl;
    if (i % n == 0)
      simul.SetOmg(0.0+1.0 - 0.2*(i/n));
  }
  simul.SetOmg(0.0+0);*/
  //simul.Dump(out_file);
  //simul.DumpForce(out_file);
  for (i = 0; i < simul.step; ++i) {
    simul.PeriodBoundary();
    simul.UpdateMNHC_VV3();
    if (!simul.ok) {
      simul.Dump(out_file);
      simul.DumpForce(out_file);
      simul.DumpEnergy(out_file);
      simul.NormalizeStat();
      simul.DumpStat(out_file);
      return 0;
    }
    simul.MakeStat();
    if (i % 100000 == 0)
      out_file << "i=" << i << std::endl;
  }
  //simul.Dump(out_file);
  simul.NormalizeStat();
  simul.DumpStat(out_file);
  t = clock() - t;
  out_file << (int)(((double)1000 * t) / CLOCKS_PER_SEC) << std::endl;
  return 0;
}
