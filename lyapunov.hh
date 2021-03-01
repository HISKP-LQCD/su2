#pragma once

#include"gaugeconfig.hh"
#include"adjointfield.hh"
#include"hamiltonian_field.hh"
#include"md_params.hh"
#include"integrator.hh"
#include<iostream>
#include<vector>
#include<list>
#include<random>
#include<iostream>
#include<fstream>
#include<cmath>
#include<complex>

using std::vector;

// this performs two parallel MD integrations
// starting from identical momenta
// and up to rounding the same gauge fields

template<class T, class URNG> void compute_lyapunov(gaugeconfig<su2> &U, 
                                                    URNG &engine, 
                                                    md_params params,
                                                    std::list<monomial<T>*> &monomial_list, 
                                                    integrator<T> &md_integ, 
                                                    std::string const &path, 
                                                    const size_t d = 12) {

  const size_t n = pow(10, d);

  adjointfield<T> momenta(U.getLs(), U.getLt());
  // generate standard normal distributed random momenta
  momenta = initnormal<URNG, T>(engine, U.getLs(), U.getLt());
  adjointfield<T> momenta2(momenta);
  
  // generate copy of U, but round to d decimal digits
  gaugeconfig<su2> U2(U.getLs(), U.getLt(), U.getBeta());
  if(d != 0) {
    for(size_t i = 0; i < U.getSize(); i++) {
      U2[i] = U[i].round(n);
    }
  }
  else {
    for(size_t i = 0; i < U.getSize(); i++) {
      U2[i] = U[i];
    }
  }
  
  hamiltonian_field<T> h(momenta, U);
  hamiltonian_field<T> h2(momenta2, U2);

  double dtau = params.gettau() / params.getnsteps();
  params.settau(dtau);
  params.setnsteps(1);

  std::ofstream os(path, std::ios::out);
  for(size_t t = 0; t < 100; t++) {
    md_integ.integrate(monomial_list, h, params, false);
    md_integ.integrate(monomial_list, h2, params, false);
    double sum = 0.;
    for(size_t i = 0; i < U.getSize(); i++) {
      sum += norm(U[i].geta() - U2[i].geta()) + norm(U[i].getb() - U2[i].getb());
    }
    os << t << " " << sum << std::endl;
  }
  return;
}
