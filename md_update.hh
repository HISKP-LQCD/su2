#pragma once

#include"gaugeconfig.hh"
#include"hamiltonian_field.hh"
#include"monomial.hh"
#include"md_update.hh"
#include<iostream>
#include<vector>
#include<list>
#include<random>

using std::vector;

class md_params {

public:
  md_params(size_t _nsteps, double _tau) : nsteps(_nsteps), tau(_tau) {}
  size_t getnsteps() const {
    return nsteps;
  }
  double gettau() const {
    return tau;
  }
  void setnsteps(const size_t _nsteps) {
    nsteps = _nsteps;
  }
  void settau(const double _tau) {
    tau = _tau;
  }

private:
  size_t nsteps;
  double tau;
};

template<class URNG> int md_update(gaugeconfig &U,
                                   URNG &engine, 
                                   md_params const &params,
                                   std::list<monomial<double>*> &monomial_list) {
  size_t Volume = U.getVolume();
  vector<double> momenta;
  momenta.resize(Volume*4*3);
  std::normal_distribution<double> normal(0., 1.);
  std::uniform_real_distribution<double> uniform(0., 1.);
  // generate standard normal distributed random momenta
  for(int i = 0; i < Volume*4*3; i++) {
    momenta[i] = normal(engine);
  }
  hamiltonian_field<double> h(momenta, U);

  // compute the initial Hamiltonian
  for (std::list<monomial<double>*>::iterator it = monomial_list.begin(); it != monomial_list.end(); it++) {
    (*it)->heatbath(h); 
  }

  gaugeconfig U_old(U);
  // perform MD evolution


  // compute the final Hamiltonian
  double delta_H = 0.;
  for (std::list<monomial<double>*>::iterator it = monomial_list.begin(); it != monomial_list.end(); it++) {
    (*it)->accept(h); 
    delta_H += (*it)->getDeltaH();
  }
  std::cout << delta_H << std::endl;

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

