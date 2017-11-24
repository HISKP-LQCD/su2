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

enum integrators { LEAPFROG = 0, LP_LEAPFROG = 1, OMF4 = 2, LP_OMF4 = 3, EULER = 4};

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
    // restore SU
    h.U->restoreSU();
  }
};

// leapfrog with rounding to allow for lower precision during MD update
template<class T> class lp_leapfrog : public integrator<T> {
public:
  lp_leapfrog() : n_prec(16) {}
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
    // restore SU
    h.U->restoreSU();
  }
private:
  size_t n_prec;
};


// OMF4 integration scheme
template<class T> class omf4 : public integrator<T> {
public:
  omf4() : rho(0.2539785108410595), theta(-0.03230286765269967), vartheta(0.08398315262876693), lambda(0.6822365335719091) {}
  void integrate(std::list<monomial<T>*> &monomial_list, hamiltonian_field<T> &h, md_params const &params) {
    adjointfield<T> deriv(h.U->getLs(), h.U->getLt());
    
    T dtau = params.gettau()/T(params.getnsteps());
    T eps[10] = {rho*dtau, lambda*dtau, 
                 theta*dtau, 0.5*(1-2.*(lambda+vartheta))*dtau, 
                 (1-2.*(theta+rho))*dtau, 0.5*(1-2.*(lambda+vartheta))*dtau,
                 theta*dtau, lambda*dtau,
                 rho*dtau, 2*vartheta*dtau};
    
    // initial half-step for the momenta
    update_momenta(monomial_list, deriv, h, 0.5*eps[9]);
    // nsteps-1 full steps
    for(size_t i = 1; i < params.getnsteps()-1; i++) {
      for(size_t j = 0; j < 5; j++) {
        update_gauge(h, eps[2*j]);
        update_momenta(monomial_list, deriv, h, eps[2*j+1]);
      }
    }
    // almost one more full step
    for(size_t j = 0; j < 4; j++) {
      update_gauge(h, eps[2*j]);
      update_momenta(monomial_list, deriv, h, eps[2*j+1]);
    }
    update_gauge(h, eps[8]);
    // final half-step in the momenta
    update_momenta(monomial_list, deriv, h, 0.5*eps[9]);
    // restore SU
    h.U->restoreSU();
  }
private:
  double rho, theta, vartheta, lambda;
};

// OMF4 integration scheme in low precision
template<class T> class lp_omf4 : public integrator<T> {
public:
  lp_omf4() : rho(0.2539785108410595), theta(-0.03230286765269967), vartheta(0.08398315262876693), lambda(0.6822365335719091), n_prec(16) {}
  lp_omf4(size_t n) : rho(0.2539785108410595), theta(-0.03230286765269967), vartheta(0.08398315262876693), lambda(0.6822365335719091), n_prec(n) {}
  void integrate(std::list<monomial<T>*> &monomial_list, hamiltonian_field<T> &h, md_params const &params) {
    adjointfield<T> deriv(h.U->getLs(), h.U->getLt());
    const size_t N = pow(10, n_prec);

    T dtau = params.gettau()/T(params.getnsteps());
    T eps[10] = {rho*dtau, lambda*dtau, 
                 theta*dtau, 0.5*(1-2.*(lambda+vartheta))*dtau, 
                 (1-2.*(theta+rho))*dtau, 0.5*(1-2.*(lambda+vartheta))*dtau,
                 theta*dtau, lambda*dtau,
                 rho*dtau, 2*vartheta*dtau};
    
    // initial half-step for the momenta
    round_and_update_momenta(monomial_list, deriv, h, 0.5*eps[9], N);
    // nsteps-1 full steps
    for(size_t i = 1; i < params.getnsteps()-1; i++) {
      for(size_t j = 0; j < 5; j++) {
        round_and_update_gauge(h, eps[2*j], N);
        round_and_update_momenta(monomial_list, deriv, h, eps[2*j+1], N);
      }
    }
    // almost one more full step
    for(size_t j = 0; j < 4; j++) {
      round_and_update_gauge(h, eps[2*j], N);
      round_and_update_momenta(monomial_list, deriv, h, eps[2*j+1], N);
    }
    round_and_update_gauge(h, eps[8], N);
    // final half-step in the momenta
    round_and_update_momenta(monomial_list, deriv, h, 0.5*eps[9], N);
    // restore SU
    h.U->restoreSU();
  }
private:
  double rho, theta, vartheta, lambda;
  size_t n_prec;
};

// semi-implicit symplectic Euler integration scheme
template<class T> class euler : public integrator<T> {
public:
  euler() {}
  void integrate(std::list<monomial<T>*> &monomial_list, hamiltonian_field<T> &h, md_params const &params) {
    adjointfield<T> deriv(h.U->getLs(), h.U->getLt());
    
    T dtau = params.gettau()/T(params.getnsteps());
    // nsteps full steps
    for(size_t i = 0; i < params.getnsteps(); i++) {
      update_momenta(monomial_list, deriv, h, dtau);
      update_gauge(h, dtau);
    }
    // restore SU
    h.U->restoreSU();
  }
};


template<class T> integrator<T>* set_integrator(const size_t integs, const size_t exponent) {
  integrator<T> * integ;
  if(static_cast<integrators>(integs) == LEAPFROG) {
    integ = new leapfrog<T>();
    std::cerr << "leapfrog" << std::endl;
  }
  else if(static_cast<integrators>(integs) == LP_LEAPFROG) {
    integ = new lp_leapfrog<T>(exponent);
    std::cerr << "lp_leapfrog" << std::endl;
  }
  else   if(static_cast<integrators>(integs) == OMF4) {
    integ = new omf4<T>();
    std::cerr << "omf4" << std::endl;
  }
  else if(static_cast<integrators>(integs) == LP_OMF4) {
    integ = new lp_omf4<T>(exponent);
    std::cerr << "lp_omf4" << std::endl;
  }
  else if(static_cast<integrators>(integs) == EULER) {
    integ = new euler<T>();
    std::cerr << "euler" << std::endl;
  }
  else {
    std::cerr << "Integrator does not match, using default" << std::endl;
    integ = new leapfrog<T>();
  }
  return integ;
}
