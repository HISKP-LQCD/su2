// Copyright (C) 2017 C. Urbach

#pragma once

#include"accum_type.hh"
#include"su2.hh"
#include"gaugeconfig.hh"
#include<vector>
#include<fstream>
#include<iomanip>

template<class Group=su2, typename lint=int> double planar_wilsonloop_dir(gaugeconfig<Group> &U, const size_t r, const size_t t, 
                                                                          const size_t mu, const size_t nu) {
  double loop = 0.;
  typedef typename accum_type<Group>::type accum;
  
  std::vector<lint> x = {0, 0, 0, 0};
  for (x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for (x[1] = 0; x[1] < U.getLx(); x[1]++) {
      for (x[2] = 0; x[2] < U.getLy(); x[2]++) {
        for (x[3] = 0; x[3] < U.getLz(); x[3]++) {
          std::vector<lint> xrun = x;
          accum L(1., 0.);
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
          loop += retrace(L);
        }
      }
    }
  }
  return loop;
}

template<class Group> double wilsonloop(gaugeconfig<Group> &U, const size_t r, const size_t t) {
  double loop = 0.;
  for(size_t mu = 1; mu < U.getndims(); mu++) {
    loop += planar_wilsonloop_dir(U, r, t, mu, 0);
  }
  return loop/U.getVolume()/double(U.getNc())/3.;
}

template<class Group> void compute_all_loops(gaugeconfig<Group> &U, std::string const &path) {
  std::ofstream os(path, std::ios::out);
  for(size_t t = 1; t < U.getLt(); t++) {
    os << t << " ";
    for(size_t r = 1; r < U.getLx(); r++) {
      double loop = wilsonloop(U, r, t);
      os << std::scientific << std::setw(15) << loop << " ";
    }
    os << std::endl;
  }
  return;
}

template<class Group> void compute_spacial_loops(gaugeconfig<Group> &U, std::string const &path) {
  std::ofstream os(path, std::ios::out);
  size_t r[2] = {2, 8};
  for(size_t t = 1; t < U.getLt(); t++) {
    os << t << " ";
    for(size_t i = 0; i < 2; i++) {
      double loop = wilsonloop(U, r[i], t);
      os << std::scientific << std::setw(15) << loop << " ";
    }
    os << std::endl;
  }
  return;
}
