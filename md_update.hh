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


template<class URNG, class T> void md_update(gaugeconfig &U,
                                             URNG &engine, 
                                             md_params &params,
                                             std::list<monomial<T>*> &monomial_list, 
                                             integrator<T> &md_integ) {
  adjointfield<T> momenta(U.getLs(), U.getLt());
  // generate standard normal distributed random momenta
  // normal distribution checked!
  momenta = initnormal<URNG, T>(engine, U.getLs(), U.getLt());

  std::uniform_real_distribution<T> uniform(0., 1.);

  hamiltonian_field<T> h(momenta, U);

  // compute the initial Hamiltonian
  for (auto it = monomial_list.begin(); it != monomial_list.end(); it++) {
    (*it)->heatbath(h); 
  }

  // keep a copy of original gauge field
  gaugeconfig U_old(U);

  // perform MD evolution
  md_integ.integrate(monomial_list, h, params);

  // compute the final Hamiltonian
  double delta_H = 0.;
  // perform acceptance part in monomials and
  // collect all the deltaH pieces from the different monomials
  for (auto it = monomial_list.begin(); it != monomial_list.end(); it++) {
    (*it)->accept(h); 
    delta_H += (*it)->getDeltaH();
  }
  params.setdeltaH(delta_H);

  // accept/reject step, if needed
  params.setaccept(true);
  if(delta_H > 0) {
    if(uniform(engine) > exp(-delta_H)) {
      params.setaccept(false);
    }
  }

  // if wanted, perform a reversibility violation test.
  if(params.getrevtest()) {
    delta_H = 0.;
    gaugeconfig U_save(U);
    h.momenta->flipsign();
    md_integ.integrate(monomial_list, h, params);

    for (auto it = monomial_list.begin(); it != monomial_list.end(); it++) {
      (*it)->accept(h); 
      delta_H += (*it)->getDeltaH();
    }
    params.setdeltadeltaH(delta_H);
    U = U_save;
  }

  // in case not accepted, restore initial gauge field
  if(!params.getaccept()) {
    U = U_old;
  }
  return;
}


