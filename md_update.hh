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


template<class URNG> void md_update(gaugeconfig &U,
                                    URNG &engine, 
                                    md_params &params,
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
  if(params.getexponent() < 1) {
    leapfrog<double> md_integ;
    md_integ.integrate(monomial_list, h, params);
  }
  else {
    lp_leapfrog<double> md_integ(params.getexponent());
    md_integ.integrate(monomial_list, h, params);
  }

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
    if(params.getexponent() < 1) {
      leapfrog<double> md_integ;
      md_integ.integrate(monomial_list, h, params);
    }
    else {
      lp_leapfrog<double> md_integ(params.getexponent());
      md_integ.integrate(monomial_list, h, params);
    }

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


