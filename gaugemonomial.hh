#pragma once
#include"monomial.hh"
#include"gaugeconfig.hh"
#include"hamiltonian_field.hh"
#include"gauge_energy.hh"

// gauge monomial
template<class T> class gaugemonomial : public monomial<T> {
public:
  gaugemonomial<T>(unsigned int _timescale) : monomial<T>::monomial(_timescale) {}
  void heatbath(hamiltonian_field<T> const &h) override {
    monomial<T>::Hold = h.U->getBeta()/N_c*gauge_energy(*(h.U));
    return;
  }
  void accept(hamiltonian_field<T> const &h) override {
    monomial<T>::Hnew = h.U->getBeta()/N_c*gauge_energy(*(h.U));
    return;
  }
  void derivative(std::vector<T> &x, hamiltonian_field<T> const &h) const override {
    return;
  }
};
