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
  leapfrog<double> md_integ;
  md_integ.integrate(monomial_list, h, params);

  // compute the final Hamiltonian
  double delta_H = 0.;
  // collect all the pieces from the different monomials
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


template<class URNG> void kramers_md_update(gaugeconfig &U,
                                            URNG &engine, 
                                            md_params &params,
                                            std::list<monomial<double>*> &monomial_list) {
  adjointfield<double> momenta(U.getLs(), U.getLt());
  // generate standard normal distributed random momenta
  // normal distribution checked!
  momenta = initnormal<URNG, double>(engine, U.getLs(), U.getLt());

  // for the accept reject step
  std::uniform_real_distribution<double> uniform(0., 1.);

  const double gamma = 2.0;

  for(size_t k = 0; k < params.getkmax(); k++) {
    // first momenta update
    
    {
      adjointfield<double> eta(U.getLs(), U.getLt());
      // generate standard normal distributed random momenta
      // normal distribution checked!
      eta = initnormal<URNG, double>(engine, U.getLs(), U.getLt());
      double epsilon = params.gettau()/params.getnsteps();
      for(size_t i = 0; i < eta.getSize(); i++) {
        momenta[i].seta(momenta[i].geta()*exp(-gamma*epsilon) + sqrt(1 - exp(-2*gamma*epsilon))*eta[i].geta());
        momenta[i].setb(momenta[i].getb()*exp(-gamma*epsilon) + sqrt(1 - exp(-2*gamma*epsilon))*eta[i].getb());
        momenta[i].setc(momenta[i].getc()*exp(-gamma*epsilon) + sqrt(1 - exp(-2*gamma*epsilon))*eta[i].getc());
      }
    }
    adjointfield<double> momenta_old(U.getLs(), U.getLt());
    momenta_old = momenta;

    // keep a copy of original gauge field and the momenta
    gaugeconfig U_old(U);
    
    hamiltonian_field<double> h(momenta, U);

    // compute the initial Hamiltonian
    for (auto it = monomial_list.begin(); it != monomial_list.end(); it++) {
      (*it)->heatbath(h); 
    }
    
    // perform MD evolution
    leapfrog<double> md_integ;
    md_integ.integrate(monomial_list, h, params);
    
    // compute the final Hamiltonian
    double delta_H = 0.;
    // collect all the pieces from the different monomials
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

    // in case not accepted, restore initial gauge field
    if(!params.getaccept()) {
      U = U_old;
      for(size_t i = 0; i < momenta.getSize(); i++) {
        momenta[i].seta(-momenta_old[i].geta());
        momenta[i].setb(-momenta_old[i].getb());
        momenta[i].setc(-momenta_old[i].getc());
      }
    }
  }
  return;
}

