#pragma once

#include"monomial.hh"
#include"adjointfield.hh"
#include"hamiltonian_field.hh"
#include<list>
#include<vector>
#include<iostream>

template<class T> void update_momenta(std::list<monomial<T>*> &monomial_list, 
                                      adjointfield<T> &deriv, hamiltonian_field<T> &h, 
                                      const double dtau) {

  zeroadjointfield(deriv);

  // compute derivatives
  // NOTE:
  // by using the intermediate field deriv we allow us to later
  // introduce more than one monomial per timescale
  for (std::list<monomial<double>*>::iterator it = monomial_list.begin(); it != monomial_list.end(); it++) {
    if((*it)->getmdactive() && ((*it)->getTimescale() == 0)) {
      (*it)->derivative(deriv, h, 1.);
    }
  }
  
  // update the momenta after the derivative has been accumulated
#pragma omp parallel for
  for(size_t i = 0; i < h.momenta->getSize(); i++) {
    (*h.momenta)[i] -= dtau * deriv[i];
  }

  return;
}

template<class T> void round_and_update_momenta(std::list<monomial<T>*> &monomial_list, 
                                                adjointfield<T> &deriv, hamiltonian_field<T> &h, 
                                                const double dtau, const size_t n) {
  
  zeroadjointfield(deriv);

  // compute derivatives
  // NOTE:
  // by using the intermediate field deriv we allow us to later
  // introduce more than one monomial per timescale
  for (std::list<monomial<double>*>::iterator it = monomial_list.begin(); it != monomial_list.end(); it++) {
    if((*it)->getmdactive() && ((*it)->getTimescale() == 0)) {
      (*it)->derivative(deriv, h, 1.);
    }
  }
  
  // update the momenta after the derivative has been accumulated
  // round before update
#pragma omp parallel for
  for(size_t i = 0; i < h.momenta->getSize(); i++) {
    (*h.momenta)[i] -= dtau * deriv[i].round(n);
  }

  return;
}
