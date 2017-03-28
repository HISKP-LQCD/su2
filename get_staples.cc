#include"su2.hh"
#include"gaugeconfig.hh"
#include"get_staples.hh"

su2 get_staples(gaugeconfig &U, vector<size_t> const x, const size_t mu) {
  su2 K(0., 0.);
  vector<size_t> x1 = x, x2 = x;
  x1[mu] += 1;
  for(int nu = 0; nu < 4; nu++) {
    if(nu != mu) {
      x2[nu] += 1;
      K += U(x1, nu) * U(x1, mu).dagger() * U(x, nu).dagger();
      x2[nu] -= 1;
    }
  }
  for(int nu = 0; nu < 4; nu++) {
    if(nu != mu) {
      x1[nu] -= 1;
      x2[nu] -= 1;
      K += U(x1, nu).dagger() * U(x2, mu).dagger() * U(x2, nu);
      x2[nu] += 1;
    }
  }
  return(K);
}
