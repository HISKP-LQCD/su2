#pragma once

#include"gaugeconfig.hh"
#include"adjointfield.hh"
#include"hamiltonian_field.hh"
#include"monomial.hh"
#include"md_params.hh"
#include"integrator.hh"
#include<iostream>
#include<vector>
#include<list>
#include<random>
#include<iostream>

using std::vector;


template<class URNG> int md_update(gaugeconfig &U,
                                   URNG &engine, 
                                   md_params const &params,
                                   std::list<monomial<double>*> &monomial_list) {
  adjointfield<double> momenta(U.getLs(), U.getLt());
  // generate standard normal distributed random momenta
  // normal distribution checked!
  momenta = initnormal<URNG, double>(engine, U.getLs(), U.getLt());

  std::uniform_real_distribution<double> uniform(0., 1.);

  hamiltonian_field<double> h(momenta, U);

  // compute the initial Hamiltonian
  for (auto it = monomial_list.begin(); it != monomial_list.end(); it++) {
    (*it)->heatbath(h); 
  }

  // keep a copy of original gauge field
  gaugeconfig U_old(U);

  // perform MD evolution
  leapfrog<double> md_integ;
  md_integ.integrate(monomial_list, h, params);

  // compute the final Hamiltonian
  double delta_H = 0.;
  for (auto it = monomial_list.begin(); it != monomial_list.end(); it++) {
    (*it)->accept(h); 
    delta_H += (*it)->getDeltaH();
    //std::cout << "monomial deltaH: " << (*it)->getDeltaH() << " ";
  }
  std::cout << "deltaH: " << delta_H << " ";// << std::endl;

  // accept/reject step, if needed
  bool accepted = true;
  if(delta_H > 0) {
    if(uniform(engine) > exp(-delta_H)) {
      accepted = false;
    }
  }
  if(!accepted) {
    U = U_old;
  }

  return accepted;
}

