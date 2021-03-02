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
  // S_g = sum_x sum_{mu<nu} beta*(1- 1/Nc*Re[Tr[U_{mu nu}]])
  // beta = 2*N_c/g_0^2
  void heatbath(hamiltonian_field<T> const &h) override {
    monomial<T>::Hold = h.U->getBeta()*(h.U->getVolume()*6 - gauge_energy(*(h.U))/N_c);
    return;
  }
  void accept(hamiltonian_field<T> const &h) override {
    monomial<T>::Hnew = h.U->getBeta()*(h.U->getVolume()*6 - gauge_energy(*(h.U))/N_c);
    return;
  }
  void derivative(adjointfield<T> &deriv, hamiltonian_field<T> const &h, const T fac = 1.) const override {
    std::vector<size_t> x = {0, 0, 0, 0};
#pragma omp parallel for
    for(size_t x0 = 0; x0 < h.U->getLt(); x0++) {
      for(size_t x1 = 0; x1 < h.U->getLx(); x1++) {
        for(size_t x2 = 0; x2 < h.U->getLy(); x2++) {
          for(size_t x3 = 0; x3 < h.U->getLz(); x3++) {
            std::vector<size_t> x = {x0, x1, x2, x3};
            for(size_t mu = 0; mu < h.U->getndims(); mu++) {
              _su2 S(0., 0.);
              get_staples(S, *h.U, x, mu);
              S = (*h.U)(x, mu) * S;
              const Complex a = S.geta(), b = S.getb();
              // the antihermitian traceless part
              // beta/N_c *(U*U^stap - (U*U^stap)^dagger)
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
