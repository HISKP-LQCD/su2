#pragma once

#include"gaugeconfig.hh"
#include<vector>

template<class T> struct hamiltonian_field {
  vector<T> * momenta;
  gaugeconfig * U;
  hamiltonian_field(vector<T> &momenta, gaugeconfig &U) :
    momenta(&momenta), U(&U) {}
};
