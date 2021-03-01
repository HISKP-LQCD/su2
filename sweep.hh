#pragma once

#include"gaugeconfig.hh"
#include"gaugeconfig.hh"
#include"random_su2.hh"
#include"get_staples.hh"
#include<random>
#include<vector>

template<class URNG, class T> double sweep(gaugeconfig<T> &U, URNG &engine,
                                           const double delta, 
                                           const size_t N_hit, const double beta) {

  std::uniform_real_distribution<double> uniform(0., 1.);

  size_t rate = 0;
  T R;
  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLs(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLs(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLs(); x[3]++) {
          for(size_t mu = 0; mu < 4; mu++) {
            T K;
            get_staples(K, U, x, mu);
            for(size_t n = 0; n < N_hit; n++) {
              random_element(R, engine, delta);
              double deltaS = beta/static_cast<double>(N_c)*
                (trace(U(x, mu) * K) - trace(U(x, mu) * R * K));
              bool accept = (deltaS < 0);
              if(!accept) accept = (uniform(engine) < exp(-deltaS));
              if(accept) {
                U(x, mu) = U(x, mu) * R;
                U(x, mu).restoreSU();
                rate += 1;
              }
            }
          }
        }
      }
    }
  }
  return( double(rate)/double(N_hit)/double(U.getSize()));
}

template<class URNG> double sweep(gaugeconfig<_u1> &U, URNG &engine,
                                  const double delta, 
                                  const size_t N_hit, const double beta) {
  
  std::uniform_real_distribution<double> uniform(0., 1.);

  size_t rate = 0;
  _u1 R;
  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLs(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLs(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLs(); x[3]++) {
          for(size_t mu = 0; mu < 4; mu++) {
            Complex K;
            get_staples(K, U, x, mu);
            for(size_t n = 0; n < N_hit; n++) {
              random_element(R, engine, delta);
              double deltaS = beta/static_cast<double>(N_c)*
                (trace(U(x, mu) * K) - trace(U(x, mu) * R * K));
              bool accept = (deltaS < 0);
              if(!accept) accept = (uniform(engine) < exp(-deltaS));
              if(accept) {
                U(x, mu) = U(x, mu) * R;
                U(x, mu).restoreSU();
                rate += 1;
              }
            }
          }
        }
      }
    }
  }
  return( double(rate)/double(N_hit)/double(U.getSize()));
}
