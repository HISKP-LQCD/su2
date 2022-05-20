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
//not tested for SU2, for U1 lattice stays the same when smearing with alpha=1
//norm is the number of staples which are used in the summation
//if spacial=true, only the staples with mu>=1 are smeared and taken into account for the staples
template<class Group> void smearlatticeape(gaugeconfig<Group> &U, double alpha, bool spatial_only=false){
  if (U.getndims() == 2 && spatial_only){
    std::cerr << "Spacial smearing is not possible in two dimensions!" << std::endl;
    return;
  }
  gaugeconfig<Group> Uold = U;
  typedef typename accum_type<Group>::type accum;
  size_t startmu = 0;
  size_t norm = 2.0*(U.getndims()-1);
  if(spatial_only){
      startmu = 1;
      norm = 2.0*(U.getndims()-2);
  }
  #ifdef _USE_OMP_
  #pragma omp parallel for
  #endif
  for(size_t x0 = 0; x0 < U.getLt(); x0++) {
    for(size_t x1 = 0; x1 < U.getLx(); x1++) {
      for(size_t x2 = 0; x2 < U.getLy(); x2++) {
        for(size_t x3 = 0; x3 < U.getLz(); x3++) {
          std::vector<size_t> x = {x0, x1, x2, x3};  
          for(size_t mu = startmu; mu < U.getndims(); mu++) {
            // K is intialized to (0,0) even if not explicitly specified
            accum K(0.0, 0.0);
            get_staples(K, Uold, x, mu, 1.0, false, spatial_only);
            Group Uprime(Uold(x, mu) * alpha + K * (1-alpha) / double(norm));
            U(x, mu) = Uprime;
            U(x, mu).restoreSU();
          }
        }
      }
    }
  }      
}

