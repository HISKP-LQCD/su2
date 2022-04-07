#pragma once

#include"gaugeconfig.hh"

//
// U^staples = \sum_{nu != mu} U_nu(x+mu) U_mu(x+nu)^dagger + U_nu(x)^dagger
//             + sum__{nu != mu} U_nu(x+mu-nu)^dagger U_mu(x-nu)^dagger U_nu(x-nu)
// 
/**
 * calculates the staples surrounding the link at position x and direction mu
 * For anisotropic lattices, the action requires weighting the temporal and spacial plaquettes differently,
 * see e.g. https://arxiv.org/abs/hep-lat/0209159. This is achieved by including this factor in the staples
 * For other cases, only the staples involving only spatial links are needed,
 * then the summation starts at nu=1 to exclude temporal links. This is done by setting spatial_only to true
 * */

template<class T, class S> void get_staples(T &K, gaugeconfig<S> &U,
                                   vector<size_t> const x,
                                   const size_t mu, const double xi=1.0, bool anisotropic=false, bool spatial_only=false) {
  size_t startnu=0;
  if(spatial_only){
      startnu=1;
  }
  vector<size_t> x1 = x, x2 = x;
  x1[mu] += 1;
  if(!anisotropic){
  for(size_t nu = startnu; nu < U.getndims(); nu++) {
    if(nu != mu) {
      x2[nu]++;
      K += U(x1, nu) * U(x2, mu).dagger() * U(x, nu).dagger();
      x2[nu]--;
    }
  }
  for(size_t nu = startnu; nu < U.getndims(); nu++) {
    if(nu != mu) {
      x1[nu]--;
      x2[nu]--;
      K += U(x1, nu).dagger() * U(x2, mu).dagger() * U(x2, nu);
      x2[nu]++;
      x1[nu]++;
    }
  }
  }
  if(anisotropic){
    double factor;
    for(size_t nu = startnu; nu < U.getndims(); nu++) {
    if(nu != mu) {
      factor=(((nu==0)||(mu==0)) ? 1.0/xi : xi);
      x2[nu]++;
      K += factor * U(x1, nu) * U(x2, mu).dagger() * U(x, nu).dagger();
      x2[nu]--;
    }
    }
    //Maybe put both loops together so only one ?: operator is needed?
    for(size_t nu = startnu; nu < U.getndims(); nu++) {
    if(nu != mu) {
      factor=(((nu==0)||(mu==0)) ? 1.0/xi : xi);
      x1[nu]--;
      x2[nu]--;
      K += factor * U(x1, nu).dagger() * U(x2, mu).dagger() * U(x2, nu);
      x2[nu]++;
      x1[nu]++;
    }
    }
  }
}


