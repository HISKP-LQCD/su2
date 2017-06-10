#pragma once

#include"monomial.hh"
#include"adjointfield.hh"
#include"hamiltonian_field.hh"
#include<list>
#include<vector>

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
      (*it)->derivative(deriv, h);
    }
  }
  
  for(size_t i = 0; i < h.U->getSize(); i++) {
    (*h.momenta)[i] -= dtau * deriv[i];
  }
  // std::vector<size_t> x = {0, 0, 0, 0};
  // for(x[0] = 0; x[0] < (h.U)->getLt(); x[0]++) {
  //   for(x[1] = 0; x[1] < (h.U)->getLs(); x[1]++) {
  //     for(x[2] = 0; x[2] < (h.U)->getLs(); x[2]++) {
  //       for(x[3] = 0; x[3] < (h.U)->getLs(); x[3]++) {
  //         for(size_t mu = 0; mu < 4; mu++) {
  //           (*h.momenta)(x, mu) -= dtau * deriv(x, mu);
  //         }
  //       }
  //     }
  //   }
  // }
  return;
}
