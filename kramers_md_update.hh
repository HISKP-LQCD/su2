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




template<class URNG, typename Float, class Group> void kramers_md_update(gaugeconfig<su2> &U,
                                                                         URNG &engine, 
                                                                         md_params &params,
                                                                         std::list<monomial<Float, Group>*> &monomial_list, 
                                                                         integrator<Float, Group> &md_integ) {
  adjointfield<Float, Group> momenta(U.getLx(), U.getLy(), U.getLz(), U.getLt(), U.getndims());
  // generate standard normal distributed random momenta
  // normal distribution checked!
  initnormal(engine, momenta);

  // for the accept reject step
  std::uniform_real_distribution<Float> uniform(0., 1.);

  const double gamma = params.getgamma();
  adjointfield<Float, Group> eta(U.getLx(), U.getLy(), U.getLz(), U.getLt(), U.getndims());
  adjointfield<Float, Group> momenta_old(U.getLx(), U.getLy(), U.getLz(), U.getLt(), U.getndims());
  gaugeconfig<su2> U_old(U.getLx(), U.getLy(), U.getLz(), U.getLt(), U.getndims(), U.getBeta());

  for(size_t k = 0; k < params.getkmax(); k++) {
    // first momenta update
    initnormal<URNG, Float>(engine, eta);
    Float epsilon = params.gettau()/Float(params.getnsteps());
    for(size_t i = 0; i < momenta.getSize(); i++) {
      momenta[i].seta(momenta[i].geta()*exp(-gamma*epsilon) + sqrt(1 - exp(-2*gamma*epsilon))*eta[i].geta());
      momenta[i].setb(momenta[i].getb()*exp(-gamma*epsilon) + sqrt(1 - exp(-2*gamma*epsilon))*eta[i].getb());
      momenta[i].setc(momenta[i].getc()*exp(-gamma*epsilon) + sqrt(1 - exp(-2*gamma*epsilon))*eta[i].getc());
    }

    momenta_old = momenta;

    // keep a copy of original gauge field and the momenta
    U_old = U;
    
    hamiltonian_field<Float, Group> h(momenta, U);

    // compute the initial Hamiltonian
    for (auto it = monomial_list.begin(); it != monomial_list.end(); it++) {
      (*it)->heatbath(h); 
    }
    
    // perform MD evolution
    if(params.getnsteps() != 1) {
      
    }
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

    // accept/reject step, if needed and/or wanted
    params.setaccept(true);
    if(params.getacceptreject()) {
      if(delta_H > 0) {
        if(uniform(engine) > exp(-delta_H)) {
          params.setaccept(false);
        }
      }
    }

    // in case not accepted, restore initial gauge field
    // and flip sign in momenta
    if(!params.getaccept()) {
      U = U_old;
      if(k < params.getkmax()-1) {
        momenta = momenta_old;
        momenta.flipsign();
      }
    }
  }
  return;
}

