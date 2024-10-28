#pragma once

#include"su2.hh"
#include "partitionings.hh"
#include "partitionings_nn.hh"
#include"random_element.hh"
#include "gaugeconfig.hh"

#include<random>


template<class T> void random_gauge_trafo(gaugeconfig<T> &U, const int seed) {
  std::mt19937 engine(seed);
 
 #ifndef partinn  
  T rU, tmp;
  #endif
  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLx(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLy(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLz(); x[3]++) {
          std::vector<size_t> xminusmu = x;
          #ifndef partinn
          random_element(rU, engine, 1);
          #endif
          for(size_t mu = 0; mu < U.getndims(); mu++) {
            #ifndef partinn
            U(x, mu) = rU * U(x, mu);
            #else
            U(x, mu) = U(x, mu);
            #endif
            xminusmu[mu] -= 1;
            #ifndef partinn
            U(xminusmu, mu) = U(xminusmu, mu) * rU.dagger();
            #else
            //This is cheating, but the only way to get this work without a group
            U(xminusmu, mu) = U(xminusmu, mu);
            #endif
            xminusmu[mu] += 1;
          }
        }
      }
    }
  }
  return;
}
