#pragma once
#include"su2.hh"
#include"monomial.hh"
#include"gaugeconfig.hh"
#include"adjointfield.hh"
#include"hamiltonian_field.hh"
#include"gauge_energy.hh"
#include"get_staples.hh"
#include<vector>
#include<complex>

// gauge monomial
template<class T> class gaugemonomial : public monomial<T> {
public:
  gaugemonomial<T>(unsigned int _timescale) : monomial<T>::monomial(_timescale) {}
  void heatbath(hamiltonian_field<T> const &h) override {
    monomial<T>::Hold = h.U->getBeta()/N_c*(h.U->getVolume()*6 - gauge_energy(*(h.U)));
    return;
  }
  void accept(hamiltonian_field<T> const &h) override {
    monomial<T>::Hnew = h.U->getBeta()/N_c*(h.U->getVolume()*6 - gauge_energy(*(h.U)));
    return;
  }
  void derivative(adjointfield<T> &deriv, hamiltonian_field<T> const &h, const T fac = 1.) const override {
    std::vector<size_t> x = {0, 0, 0, 0};
    for(x[0] = 0; x[0] < h.U->getLt(); x[0]++) {
      for(x[1] = 0; x[1] < h.U->getLs(); x[1]++) {
        for(x[2] = 0; x[2] < h.U->getLs(); x[2]++) {
          for(x[3] = 0; x[3] < h.U->getLs(); x[3]++) {
            for(size_t mu = 0; mu < 4; mu++) {
              _su2 S = (*h.U)(x, mu) * get_staples(*h.U, x, mu);
              const Complex a = S.geta(), b = S.getb();
              deriv(x, mu) += fac*h.U->getBeta()/double(N_c) *
                adjoint<double>(2.*std::imag(b), 2.*std::real(b), 2.*std::imag(a));
            }
          }
        }
      }
    }
    return;
  }
};
