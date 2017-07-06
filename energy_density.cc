#include"su2.hh"
#include"gaugeconfig.hh"
#include"energy_density.hh"
#include<complex>

// symmetric definition of the energy density
// using the clover leaf
//
//         
//  <--   <--
// |  ^  |  ^
// |  |  |  | mu
// -->   -->
//
//  <--   <--
// |  ^  |  ^
// |  |  |  |
// -->   -->
//        nu
//
// checked for gauge invariance

double energy_density(gaugeconfig &U) {
  double res = 0.;

  std::vector<size_t> x = {0, 0, 0, 0};
  for(x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for(x[1] = 0; x[1] < U.getLs(); x[1]++) {
      for(x[2] = 0; x[2] < U.getLs(); x[2]++) {
        for(x[3] = 0; x[3] < U.getLs(); x[3]++) {
          std::vector<size_t> x1 = x;
          std::vector<size_t> x2 = x;
          std::vector<size_t> x3 = x;
          for(size_t mu = 0; mu < 3; mu++) {
            for(size_t nu = mu+1; nu < 4; nu++) {
              x1[mu] += 1;
              x2[nu] += 1;
              su2 leaf = U(x, mu) * U(x1, nu) *
                U(x2, mu).dagger()*U(x, nu).dagger();
              x1[mu] -= 1;
              x2[nu] -= 1;

              x1[mu] -= 1;
              x1[nu] += 1;
              x2[mu] -= 1;
              leaf += U(x, nu) * U(x1, mu).dagger() *
                U(x2, nu).dagger()*U(x2, mu);
              x1[mu] += 1;
              x1[nu] -= 1;
              x2[mu] += 1;

              x1[mu] -= 1;
              x2[mu] -= 1;
              x2[nu] -= 1;
              x3[nu] -= 1;
              leaf += U(x1, mu).dagger() * U(x2, nu).dagger() *
                U(x2, mu)*U(x3, nu);
              x1[mu] += 1;
              x2[mu] += 1;
              x2[nu] += 1;
              x3[nu] += 1;
              
              x1[nu] -= 1;
              x2[nu] -= 1;
              x2[mu] += 1;
              leaf += U(x1, nu).dagger() * U(x1, mu) *
                U(x2, nu)*U(x, mu).dagger();
              x1[nu] += 1;
              x2[nu] += 1;
              x2[mu] -= 1;

              // traceless and anti-hermitian
              su2 one(0.5*(leaf.geta()-std::conj(leaf.geta())), 
                      leaf.getb());

              res += trace(one*one);
            }
          }
        }
      }
    }
  }
  return(res/U.getVolume());
}
