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

enum integrators { LEAPFROG = 0, LP_LEAPFROG = 1, OMF4 = 2, LP_OMF4 = 3, EULER = 4, RUTH = 5, OMF2 = 6};

// virtual integrator class
template<typename Float, class Group> class integrator{
public:
  integrator() {}
  virtual void integrate(std::list<monomial<Float, Group>*> &monomial_list, hamiltonian_field<Float, Group> &h, 
                         md_params const &params, const bool restore = true) = 0;
};

// leapfrog integration scheme
template<typename Float, class Group> class leapfrog : public integrator<Float, Group> {
public:
  leapfrog() {}
  void integrate(std::list<monomial<Float, Group>*> &monomial_list, hamiltonian_field<Float, Group> &h, 
                 md_params const &params, const bool restore=true) {

    adjointfield<typename adjoint_type<Float, Group>::type> deriv(h.U->getLx(), h.U->getLy(), h.U->getLz(), h.U->getLt(), h.U->getndims());
    
    Float dtau = params.gettau()/Float(params.getnsteps());
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
    if(restore) h.U->restoreSU();
  }
};

// leapfrog with rounding to allow for lower precision during MD update
template<typename Float, class Group> class lp_leapfrog : public integrator<Float, Group> {
public:
  lp_leapfrog() : n_prec(16) {}
  lp_leapfrog(size_t n) : n_prec(n) {}
  void integrate(std::list<monomial<Float, Group>*> &monomial_list, hamiltonian_field<Float, Group> &h, 
                 md_params const &params, const bool restore = true) {
    adjointfield<typename adjoint_type<Float, Group>::type> deriv(h.U->getLx(), h.U->getLy(), h.U->getLz(), h.U->getLt(), h.U->getndims());
    const size_t N = pow(10, n_prec);

    Float dtau = params.gettau()/Float(params.getnsteps());
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
    if(restore) h.U->restoreSU();
  }
private:
  size_t n_prec;
};

// OMF2 integration scheme (also called 2MN, second order minimal norm)
template<typename Float, class Group> class omf2 : public integrator<Float, Group> {
public:
  omf2() : lambda(0.1938) {}
  void integrate(std::list<monomial<Float, Group>*> &monomial_list, hamiltonian_field<Float, Group> &h, 
                 md_params const &params, const bool restore = true) {

    adjointfield<typename adjoint_type<Float, Group>::type> deriv(h.U->getLx(), h.U->getLy(), h.U->getLz(), h.U->getLt(), h.U->getndims());
    Float dtau = params.gettau()/Float(params.getnsteps());
    Float oneminus2lambda = (1.-2.*lambda);

    
    // initial step for the momenta
    update_momenta(monomial_list, deriv, h, lambda*dtau);
    // nsteps-1 full steps
    for(size_t i = 0; i < params.getnsteps()-1; i++) {
      update_gauge(h, 0.5*dtau);
      update_momenta(monomial_list, deriv, h, oneminus2lambda*dtau);
      update_gauge(h, 0.5*dtau);
      // double step for the momenta to avoid double computation
      update_momenta(monomial_list, deriv, h, 2*lambda*dtau);
    }
    // almost one more full step
    update_gauge(h, 0.5*dtau);
    update_momenta(monomial_list, deriv, h, oneminus2lambda*dtau);
    update_gauge(h, 0.5*dtau);
    // final step in the momenta
    update_momenta(monomial_list, deriv, h, lambda*dtau);
    // restore SU
    if(restore) h.U->restoreSU();
  }
private:
  double lambda;
};


// OMF4 integration scheme
template<typename Float, class Group> class omf4 : public integrator<Float, Group> {
public:
  omf4() : rho(0.2539785108410595), theta(-0.03230286765269967), vartheta(0.08398315262876693), lambda(0.6822365335719091) {}
  void integrate(std::list<monomial<Float, Group>*> &monomial_list, hamiltonian_field<Float, Group> &h, 
                 md_params const &params, const bool restore = true) {

    adjointfield<typename adjoint_type<Float, Group>::type> deriv(h.U->getLx(), h.U->getLy(), h.U->getLz(), h.U->getLt(), h.U->getndims());
    Float dtau = params.gettau()/Float(params.getnsteps());
    Float eps[10] = {rho*dtau, lambda*dtau, 
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
    if(restore) h.U->restoreSU();
  }
private:
  double rho, theta, vartheta, lambda;
};

// OMF4 integration scheme in low precision
template<typename Float, class Group> class lp_omf4 : public integrator<Float, Group> {
public:
  lp_omf4() : rho(0.2539785108410595), theta(-0.03230286765269967), vartheta(0.08398315262876693), lambda(0.6822365335719091), n_prec(16) {}
  lp_omf4(size_t n) : rho(0.2539785108410595), theta(-0.03230286765269967), vartheta(0.08398315262876693), lambda(0.6822365335719091), n_prec(n) {}
  void integrate(std::list<monomial<Float, Group>*> &monomial_list, hamiltonian_field<Float, Group> &h, 
                 md_params const &params, const bool restore = true) {

    adjointfield<typename adjoint_type<Float, Group>::type> deriv(h.U->getLx(), h.U->getLy(), h.U->getLz(), h.U->getLt(), h.U->getndims());
    const size_t N = pow(10, n_prec);

    Float dtau = params.gettau()/Float(params.getnsteps());
    Float eps[10] = {rho*dtau, lambda*dtau, 
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
    if(restore) h.U->restoreSU();
  }
private:
  double rho, theta, vartheta, lambda;
  size_t n_prec;
};

// semi-implicit symplectic Euler integration scheme
template<typename Float, class Group> class euler : public integrator<Float, Group> {
public:
  euler() {}
  void integrate(std::list<monomial<Float, Group>*> &monomial_list, hamiltonian_field<Float, Group> &h, 
                 md_params const &params, const bool restore = true) {

    adjointfield<typename adjoint_type<Float, Group>::type> deriv(h.U->getLx(), h.U->getLy(), h.U->getLz(), h.U->getLt(), h.U->getndims());
    Float dtau = params.gettau()/Float(params.getnsteps());
    // nsteps full steps
    for(size_t i = 0; i < params.getnsteps(); i++) {
      update_momenta(monomial_list, deriv, h, dtau);
      update_gauge(h, dtau);
    }
    // restore SU
    if(restore) h.U->restoreSU();
  }
};

// third order symplectic, but non-reversible integrator (Ruth)
template<typename Float, class Group> class ruth : public integrator<Float, Group> {
public:
  ruth() {}
  void integrate(std::list<monomial<Float, Group>*> &monomial_list, hamiltonian_field<Float, Group> &h, 
                 md_params const &params, const bool restore = true) {

    adjointfield<typename adjoint_type<Float, Group>::type> deriv(h.U->getLx(), h.U->getLy(), h.U->getLz(), h.U->getLt(), h.U->getndims());
    
    Float dtau = params.gettau()/Float(params.getnsteps());
    // nsteps full steps
    for(size_t i = 0; i < params.getnsteps(); i++) {
      update_momenta(monomial_list, deriv, h, dtau);
      update_gauge(h, -1./24.*dtau);
      update_momenta(monomial_list, deriv, h, -2./3.*dtau);
      update_gauge(h, 3./4.*dtau);
      update_momenta(monomial_list, deriv, h, 2./3.*dtau);
      update_gauge(h, 7./24.*dtau);
    }
    // restore SU
    if(restore) h.U->restoreSU();
  }
};


template<typename Float, class Group> integrator<Float, Group>* set_integrator(const size_t integs, const size_t exponent) {
  integrator<Float, Group> * integ;
  if(static_cast<integrators>(integs) == LEAPFROG) {
    integ = new leapfrog<Float, Group>();
    std::cerr << "leapfrog" << std::endl;
  }
  else if(static_cast<integrators>(integs) == LP_LEAPFROG) {
    integ = new lp_leapfrog<Float, Group>(exponent);
    std::cerr << "lp_leapfrog" << std::endl;
  }
  else if(static_cast<integrators>(integs) == OMF4) {
    integ = new omf4<Float, Group>();
    std::cerr << "omf4" << std::endl;
  }
  else if(static_cast<integrators>(integs) == LP_OMF4) {
    integ = new lp_omf4<Float, Group>(exponent);
    std::cerr << "lp_omf4" << std::endl;
  }
  else if(static_cast<integrators>(integs) == EULER) {
    integ = new euler<Float, Group>();
    std::cerr << "euler" << std::endl;
  }
  else if(static_cast<integrators>(integs) == RUTH) {
    integ = new ruth<Float, Group>();
    std::cerr << "ruth" << std::endl;
  }
  else if(static_cast<integrators>(integs) == OMF2) {
    integ = new omf2<Float, Group>();
    std::cerr << "omf2" << std::endl;
  }
  else {
    std::cerr << "Integrator does not match, using default" << std::endl;
    integ = new leapfrog<Float, Group>();
  }
  return integ;
}
