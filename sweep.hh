#pragma once

#include"su2.hh"
#include"gaugeconfig.hh"
#include"random_su2.hh"
#include"get_staples.hh"
#include<random>
#include<vector>
#include<iostream>

template<class URNG> double sweep(gaugeconfig<su2> &U, URNG &engine,
                                           const double delta, 
                                           const size_t N_hit, const double beta) {

  std::uniform_real_distribution<double> uniform(0., 1.);

  size_t rate = 0;
  su2 R(0., 0.);
  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLx(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLy(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLz(); x[3]++) {
          for(size_t mu = 0; mu < U.getndims(); mu++) {
            su2 K = get_staples(U, x, mu);
            for(size_t n = 0; n < N_hit; n++) {
              random_su2(R, engine, delta);
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

template<class URNG> double sweep(gaugeconfig<Gsu2> &U, URNG &engine,
                                  const size_t m, const size_t delta,
                                  const size_t N_hit, const double beta) {

  std::uniform_real_distribution<double> uniform(0., 1.);
  std::uniform_int_distribution<int> uniindex(0, 3);
  std::uniform_int_distribution<int> unij(0, delta-1);
  std::uniform_int_distribution<int> unisign(0, 1);
  
  size_t rate = 0;
  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLx(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLy(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLz(); x[3]++) {
          for(size_t mu = 0; mu < U.getndims(); mu++) {
            su2 K = get_staples(U, x, mu);
            for(size_t n = 0; n < N_hit; n++) {
              Gsu2 R = U(x, mu);
              // get a pair
              size_t p[2];
              p[0] = uniindex(engine);
              p[1] = p[0];
              while(p[0] == p[1]) {
                p[1] = uniindex(engine);
              }
              // get the corresponding j-values
              size_t j[2];
              j[0] = U(x,mu).getj(p[0]);
              j[1] = U(x,mu).getj(p[1]);
              if((j[0] == 0 && j[1] == 0) || (j[0] == m && j[1] == m)) {
                // only flip signs
                R.sets(p[0], (2*unisign(engine) - 1));
                R.sets(p[1], (2*unisign(engine) - 1));
              }
              else {
                bool done = false;
                while(!done) {
                  // get a shift
                  int shift = (2*unisign(engine) - 1)*(unij(engine) + 1);
                  if(j[0] + shift < 0) { // shift < 0
                    done = true;
                    size_t _j = - shift - j[0];
                    j[1] += j[0] - _j;
                    j[0] = _j;
                    // flipt sign in s[p[0]]
                    R.sets(p[0], -R.gets(p[0]));
                  }
                  else if(j[1] - shift < 0) { // shift > 0
                    done = true;
                    size_t _j = shift - j[1];
                    j[0] += j[1] - _j;
                    j[1] = _j;
                    // flip sign in s[p[1]]
                    R.sets(p[1], -R.gets(p[1]));
                  }
                  else if((j[0] + shift <= m) && (j[0] + shift >= 0) &&
                          (j[1] - shift <= m) && (j[1] - shift >= 0)) {
                    done = true;
                    j[0] += shift;
                    j[1] -= shift;
                  }
                }
                R.setjpair(p[0], p[1], j[0], j[1]);
                R.restoreSU();
              }
              double deltaS = beta/static_cast<double>(N_c)*
                (trace(U(x, mu) * K) - trace(R * K));
              bool accept = (deltaS < 0);
              if(!accept) accept = (uniform(engine) < exp(-deltaS));
              if(accept) {
                U(x, mu) = R;
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
