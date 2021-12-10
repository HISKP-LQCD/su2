#pragma once

#include"su2.hh"
#include"u1.hh"
#include"random_element.hh"
#include"gaugeconfig.hh"
#include"get_staples.hh"
#include<random>
#include<vector>
#include<cmath>
#include<fstream>
#include<complex>
#include<iostream>
#include<cassert>
#ifdef _USE_OMP_
#  include<omp.h>
#endif


//prescription for smearing taken from https://arxiv.org/abs/hep-lat/0209159
//not tested for SU2, for U1 lattice stays the same whEn smearing with alpha=1
template<class Group> void smearlatticeape(gaugeconfig<Group> &U, double alpha){
  gaugeconfig<Group> Uold = U;
  typedef typename accum_type<Group>::type accum;
  #ifdef _USE_OMP_
  #pragma omp parallel for
  #endif
  for(size_t x0 = 0; x0 < U.getLt(); x0++) {
    for(size_t x1 = 0; x1 < U.getLx(); x1++) {
      for(size_t x2 = 0; x2 < U.getLy(); x2++) {
        for(size_t x3 = 0; x3 < U.getLz(); x3++) {
          std::vector<size_t> x = {x0, x1, x2, x3};  
          for(size_t mu = 0; mu < U.getndims(); mu++) {
            accum K;
            get_staples(K, Uold, x, mu);
            K = K*(1-alpha)/2.0;
            Group Uprime(Uold(x, mu)*alpha+K);
            U(x, mu) = Uprime;
            U(x, mu).restoreSU();
          }
        }
      }
    }
  }      
}
