// Copyright (C) 2017 C. Urbach <urbach@hiskp.uni-bonn.de>

#include"wilsonloop.hh"
#include"su2.hh"
#include"gaugeconfig.hh"
#include<vector>
#include<fstream>
#include<iomanip>

double planar_wilsonloop_dir(gaugeconfig &U, const size_t r, const size_t t, 
                             const size_t mu, const size_t nu) {
  double loop = 0.;

  std::vector<size_t> x = {0, 0, 0, 0};
  for (x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for (x[1] = 0; x[1] < U.getLs(); x[1]++) {
      for (x[2] = 0; x[2] < U.getLs(); x[2]++) {
        for (x[3] = 0; x[3] < U.getLs(); x[3]++) {
          std::vector<size_t> xrun = x;
          su2 L(1., 0.);
          for (size_t _t = 0; _t < t; _t++) {
            L *= U(xrun, nu);
            xrun[nu] += 1;
          }
          for (size_t s = 0; s < r; s++) {
            L *= U(xrun, mu);
            xrun[mu] += 1;
          }
          for (size_t _t = 0; _t < t; _t++) {
            xrun[nu] -= 1;
            L *= U(xrun, nu).dagger();
          }
          for (size_t s = 0; s < r; s++) {
            xrun[mu] -= 1;
            L *= U(xrun, mu).dagger();
          }
          loop += trace(L);
        }
      }
    }
  }
  return loop;
}

double wilsonloop(gaugeconfig &U, const size_t r, const size_t t) {
  double loop = 0.;
  for(size_t mu = 1; mu < 4; mu++) {
    loop += planar_wilsonloop_dir(U, r, t, mu, 0);
  }
  return loop/U.getVolume()/N_c/3.;
}

void compute_all_loops(gaugeconfig &U, std::string const &path) {
  std::ofstream os(path, std::ios::out);
  for(size_t t = 1; t < U.getLt(); t++) {
    os << t << " ";
    for(size_t r = 1; r < U.getLs(); r++) {
      double loop = wilsonloop(U, r, t);
      os << std::scientific << std::setw(15) << loop << " ";
    }
    os << std::endl;
  }
  return;
}
