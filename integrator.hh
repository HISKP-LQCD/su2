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
#include<cmath>

// virtual integrator class
template<class T> class integrator{
public:
  integrator() {}
  virtual void integrate(std::list<monomial<T>*> &monomial_list, hamiltonian_field<T> &h, md_params const &params) = 0;
};

// leapfrog integration scheme
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
};

// leapfrog with rounding to allow for lower precision during MD update
template<class T> class lp_leapfrog : public integrator<T> {
public:
  lp_leapfrog() : n_prec(1000000000) {}
  lp_leapfrog(size_t n) : n_prec(n) {}
  void integrate(std::list<monomial<T>*> &monomial_list, hamiltonian_field<T> &h, 
                          md_params const &params) {
    adjointfield<T> deriv(h.U->getLs(), h.U->getLt());
    const size_t N = pow(10, n_prec);

    T dtau = params.gettau()/T(params.getnsteps());
    // initial half-step for the  momenta
    round_and_update_momenta(monomial_list, deriv, h, dtau/2., N);
    // first full step for gauge
    round_and_update_gauge(h, dtau, N);
    // nsteps-1 full steps
    for(size_t i = 0; i < params.getnsteps()-1; i++) {
      round_and_update_momenta(monomial_list, deriv, h, dtau, N);
      round_and_update_gauge(h, dtau, N);
    }

    // final half-step for the momenta
    round_and_update_momenta(monomial_list, deriv, h, dtau/2., N);
  }
private:
  size_t n_prec;
};
