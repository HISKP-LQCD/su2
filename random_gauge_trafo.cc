#include"su2.hh"
#include"random_su2.hh"
#include"random_gauge_trafo.hh"
#include<random>


void random_gauge_trafo(gaugeconfig &U, const int seed) {
  std::mt19937 engine(seed);
  
  su2 rU(0., 0.), tmp(0., 0.);
  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLs(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLs(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLs(); x[3]++) {
          std::vector<size_t> xminusmu = x;
          random_su2(rU, engine, 1);
          for(size_t mu = 0; mu < 4; mu++) {
            U(x, mu) = rU * U(x, mu);

            xminusmu[mu] -= 1;
            U(xminusmu, mu) = U(xminusmu, mu) * rU.dagger();
            xminusmu[mu] += 1;
          }
        }
      }
    }
  }
  return;
}
