#pragma once

#include"gaugeconfig.hh"

//
// U^staples = \sum_{nu != mu} U_nu(x+mu) U_mu(x+nu)^dagger + U_nu(x)^dagger
//             + sum__{nu != mu} U_nu(x+mu-nu)^dagger U_mu(x-nu)^dagger U_nu(x-nu)
// 

template<class T, class S> void get_staples(T &K, gaugeconfig<S> &U,
                                   vector<size_t> const x,
                                   const size_t mu) {
  vector<size_t> x1 = x, x2 = x;
  x1[mu] += 1;
  for(size_t nu = 0; nu < U.getndims(); nu++) {
    if(nu != mu) {
      x2[nu]++;
      K += U(x1, nu) * U(x2, mu).dagger() * U(x, nu).dagger();
      x2[nu]--;
    }
  }
  for(size_t nu = 0; nu < U.getndims(); nu++) {
    if(nu != mu) {
      x1[nu]--;
      x2[nu]--;
      K += U(x1, nu).dagger() * U(x2, mu).dagger() * U(x2, nu);
      x2[nu]++;
      x1[nu]++;
    }
  }
}

//stores c1*U1+complex(c2*U2) in U1
//variable K used for accumulation is of type Complex for u1 (and su2 for su2?)
inline void scalarmultiplyadd(Complex &U1, _u1 U2, double c1, double c2) {
    U1=Complex(c1,0.)*U1+Complex(c2,0.)*Complex(U2);
}

//stores c1*U1+c2*U2 in U1
void scalarmultiplyadd(_su2 &U1, _su2 U2, double c1, double c2) {
    _su2 copyU1, copyU2;
    copyU1.set(c1*U1.geta(), c1*U1.getb());
    copyU2.set(c2*U2.geta(), c2*U2.getb());
    U1=copyU1+copyU2;
}

template<class T, class S> void get_staples_anisotrope(T &K, gaugeconfig<S> &U,
                                   vector<size_t> const x,
                                   const size_t mu, const double xi=1.0) {
  vector<size_t> x1 = x, x2 = x;
  double factor=0;
  x1[mu] += 1;
  for(size_t nu = 0; nu < U.getndims(); nu++) {
    if(nu != mu) {
      factor=((nu==0)||(mu==0)?1.0/xi:xi);
      x2[nu]++;
      scalarmultiplyadd(K, U(x1, nu) * U(x2, mu).dagger() * U(x, nu).dagger(), 1.0, factor);
      x2[nu]--;
    }
  }
  for(size_t nu = 0; nu < U.getndims(); nu++) {
    if(nu != mu) {
      factor=((nu==0)||(mu==0)?1.0/xi:xi);
      x1[nu]--;
      x2[nu]--;
      scalarmultiplyadd(K, U(x1, nu).dagger() * U(x2, mu).dagger() * U(x2, nu), 1.0, factor);
      x2[nu]++;
      x1[nu]++;
    }
  }
}

