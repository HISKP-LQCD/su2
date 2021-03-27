#pragma once

#include"gaugeconfig.hh"
#include"adjointfield.hh"
#include<vector>

template<typename Float, class Group> struct hamiltonian_field {
  adjointfield<Float, Group> * momenta;
  gaugeconfig<Group> * U;
  hamiltonian_field(adjointfield<Float, Group> &momenta, gaugeconfig<Group> &U) :
    momenta(&momenta), U(&U) {}
};
