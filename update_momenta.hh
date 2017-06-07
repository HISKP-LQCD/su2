#pragma once

#include"monomial.hh"
#include"hamiltonian_field.hh"
#include<list>
#include<vector>

template<class T> void update_momenta(std::list<monomial<T>*> &monomial_list, 
                                      std::vector<T> &x, hamiltonian_field<T> &h, 
                                      const double dtau) {
  // compute derivatives
  for (std::list<monomial<double>*>::iterator it = monomial_list.begin(); it != monomial_list.end(); it++) {
    if((*it)->getmdactive() && ((*it)->getTimescale() == 0)) {
      (*it)->derivative(x, h);
    }
  }
  
  

  // set to zero
  for (auto it = x.begin(); it != x.end(); it++) {
    *it = 0.;
  }
  return;
}
