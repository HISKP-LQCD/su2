#pragma once

#include"gaugeconfig.hh"
#include"adjointfield.hh"
#include<vector>

template<class T> struct hamiltonian_field {
  adjointfield<T> * momenta;
  gaugeconfig * U;
  hamiltonian_field(adjointfield<T> &momenta, gaugeconfig &U) :
    momenta(&momenta), U(&U) {}
};
