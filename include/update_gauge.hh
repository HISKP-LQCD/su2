
//#ifndef Genz
#pragma once

#include"hamiltonian_field.hh"
#include"su2.hh"
#include "partitionings.hh"
#include"exp_gauge.hh"
#include<complex>

template<typename Float, class Group> void update_gauge(hamiltonian_field<Float, Group> &h, const Float dtau) {
  
  // update the gauge field
#pragma omp parallel for
  for(size_t i = 0; i < h.U->getSize(); i++) {
    (*h.U)[i] = exp(dtau * (*h.momenta)[i]) * (*h.U)[i];
  }
  return;
}

template<typename Float, class Group> void round_and_update_gauge(hamiltonian_field<Float, Group> &h, const Float dtau, const size_t n) {
  
  // update the gauge field
  // round before update
#pragma omp parallel for
  for(size_t i = 0; i < h.U->getSize(); i++) {
    (*h.U)[i] = exp(dtau * (*h.momenta)[i]) * (*h.U)[i].round(n);
  }
  return;
}
//#endif