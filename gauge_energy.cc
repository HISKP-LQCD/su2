#include"su2.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"


double gauge_energy(gaugeconfig &config) {
  double res = 0.;
  for(size_t t = 0; t < config.getLt(); t++) {
    for(size_t x = 0; x < config.getLs(); x++) {
      for(size_t y = 0; y < config.getLs(); y++) {
        for(size_t z = 0; z < config.getLs(); z++) {
          for(size_t mu = 0; mu < 3; mu++) {
            for(size_t nu = mu+1; nu < 4; nu++) {
              std::vector<size_t> x1 = {t, x, y, z};
              std::vector<size_t> x2 = {t, x, y, z};
              x1[mu] += 1;
              x2[nu] += 1;
              res += trace(config(t, x, y, z, mu) * config(x1, nu) *
                           config(x2, mu).dagger()*config(t, x, y, z, nu));
            }
          }
        }
      }
    }
  }
  return(res);
}
