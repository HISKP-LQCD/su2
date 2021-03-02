#pragma once

#include"gaugeconfig.hh"
#include"adjointfield.hh"
#include<vector>

template<class T> struct hamiltonian_field {
  adjointfield<T> * momenta;
  gaugeconfig<su2> * U;
  hamiltonian_field(adjointfield<T> &momenta, gaugeconfig<su2> &U) :
    momenta(&momenta), U(&U) {}
};
