#pragma once

#include"hamiltonian_field.hh"
#include"monomial.hh"
#include"md_params.hh"
#include"update_gauge.hh"
#include"update_momenta.hh"
#include"gaugeconfig.hh"
#include<vector>
#include<list>
#include<iostream>

template<class T> class integrator{
public:
  integrator() {}
  virtual void integrate(std::list<monomial<T>*> &monomial_list, hamiltonian_field<T> &h, md_params const &params) = 0;
  virtual void low_prec_integrate(std::list<monomial<T>*> &monomial_list, hamiltonian_field<T> &h, md_params const &params) = 0;
};

template<class T> class leapfrog : public integrator<T> {
public:
  leapfrog() {}
  void integrate(std::list<monomial<T>*> &monomial_list, hamiltonian_field<T> &h, md_params const &params) {
    adjointfield<T> deriv(h.U->getLs(), h.U->getLt());
    
    T dtau = params.gettau()/T(params.getnsteps());
    // initial half-step for the  momenta
    update_momenta(monomial_list, deriv, h, dtau/2.);
    // first full step for gauge
    update_gauge(h, dtau);
    // nsteps-1 full steps
    for(size_t i = 0; i < params.getnsteps()-1; i++) {
      update_momenta(monomial_list, deriv, h, dtau);
      update_gauge(h, dtau);
    }

    // final half-step for the momenta
    update_momenta(monomial_list, deriv, h, dtau/2.);
  }
  void low_prec_integrate(std::list<monomial<T>*> &monomial_list, hamiltonian_field<T> &h, 
                          md_params const &params) {
    adjointfield<T> deriv(h.U->getLs(), h.U->getLt());
    
    const size_t n = 1000000000;
    T dtau = params.gettau()/T(params.getnsteps());
    // initial half-step for the  momenta
    round_and_update_momenta(monomial_list, deriv, h, dtau/2., n);
    // first full step for gauge
    round_and_update_gauge(h, dtau, n);
    // nsteps-1 full steps
    for(size_t i = 0; i < params.getnsteps()-1; i++) {
      round_and_update_momenta(monomial_list, deriv, h, dtau, n);
      round_and_update_gauge(h, dtau, n);
    }

    // final half-step for the momenta
    round_and_update_momenta(monomial_list, deriv, h, dtau/2., n);
  }

};
