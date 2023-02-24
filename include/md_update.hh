#pragma once

#include "adjointfield.hh"
#include "gaugeconfig.hh"
#include "hamiltonian_field.hh"
#include "integrator.hh"
#include "md_params.hh"
#include "monomial.hh"
#include <iostream>
#include <list>
#include <random>
#include <vector>

using std::vector;

template <class URNG, typename Float, class Group>
void md_update(gaugeconfig<Group> &U,
               URNG &engine,
               md_params &params,
               std::list<monomial<Float, Group> *> &monomial_list,
               integrator<Float, Group> &md_integ) {
  adjointfield<Float, Group> momenta(U.getLx(), U.getLy(), U.getLz(), U.getLt(),
                                     U.getndims());
  // generate standard normal distributed random momenta
  // normal distribution checked!
  initnormal(engine, momenta);

  std::uniform_real_distribution<Float> uniform(0., 1.);

  hamiltonian_field<Float, Group> h(momenta, U);

  // compute the initial Hamiltonian
  for (auto it = monomial_list.begin(); it != monomial_list.end(); it++) {
    (*it)->heatbath(h);
  }

  // keep a copy of original gauge field
  gaugeconfig<Group> U_old(U);

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
  if (delta_H > 0) {
    if (uniform(engine) > exp(-delta_H)) {
      params.setaccept(false);
    }
  }

  // if wanted, perform a reversibility violation test.
  if (params.getrevtest()) {
    delta_H = 0.;
    gaugeconfig<Group> U_save(U);
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
  if (!params.getaccept()) {
    U = U_old;
  }
  return;
}
