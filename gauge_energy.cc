#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"su2.hh"


double gauge_energy(gaugeconfig &U) {
  double res = 0.;

  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLs(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLs(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLs(); x[3]++) {
          std::vector<size_t> xplusmu = x;
          std::vector<size_t> xplusnu = x;
          for(size_t mu = 0; mu < 3; mu++) {
            for(size_t nu = mu+1; nu < 4; nu++) {
              xplusmu[mu] += 1;
              xplusnu[nu] += 1;
              res += trace(U(x, mu) * U(xplusmu, nu) *
                           U(xplusnu, mu).dagger()*U(x, nu).dagger());
              
              xplusmu[mu] -= 1;
              xplusnu[nu] -= 1;
            }
          }
        }
      }
    }
  }
  return(res);
}
