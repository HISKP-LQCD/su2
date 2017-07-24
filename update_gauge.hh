#pragma once

#include"hamiltonian_field.hh"
#include"su2.hh"
#include"expsu2.hh"
#include<complex>

template<class T> void update_gauge(hamiltonian_field<T> &h, const T dtau) {
  
  // update the gauge field
#pragma omp parallel for
  for(size_t i = 0; i < h.U->getSize(); i++) {
    (*h.U)[i] = exp(dtau * (*h.momenta)[i]) * (*h.U)[i];
  }
  return;
}

template<class T> void round_and_update_gauge(hamiltonian_field<T> &h, const T dtau, const size_t n) {
  
  // update the gauge field
  // round before update
#pragma omp parallel for
  for(size_t i = 0; i < h.U->getSize(); i++) {
    (*h.U)[i] = exp(dtau * (*h.momenta)[i]) * (*h.U)[i].round(n);
  }
  return;
}
