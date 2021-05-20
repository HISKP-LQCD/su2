#pragma once
#include"su2.hh"
#include"u1.hh"
#include"monomial.hh"
#include"gaugeconfig.hh"
#include"adjointfield.hh"
#include"hamiltonian_field.hh"
#include"gauge_energy.hh"
#include"get_staples.hh"
#include<vector>
#include<complex>

// gauge monomial
template<typename Float, class Group, typename lint=int> class gaugemonomial : public monomial<Float, Group> {
public:
  gaugemonomial<Float, Group, lint>(unsigned int _timescale) : monomial<Float, Group>::monomial(_timescale) {}
  // S_g = sum_x sum_{mu<nu} beta*(1- 1/Nc*Re[Tr[U_{mu nu}]])
  // beta = 2*N_c/g_0^2
  void heatbath(hamiltonian_field<Float, Group> const &h) override {
    monomial<Float, Group>::Hold = h.U->getBeta()*(h.U->getVolume()*6 - gauge_energy(*(h.U))/double(h.U->getNc()));
    return;
  }
  void accept(hamiltonian_field<Float, Group> const &h) override {
    monomial<Float, Group>::Hnew = h.U->getBeta()*(h.U->getVolume()*6 - gauge_energy(*(h.U))/double(h.U->getNc()));
    return;
  }
  void derivative(adjointfield<Float, Group> &deriv, hamiltonian_field<Float, Group> const &h, const Float fac = 1.) const override {
    std::vector<lint> x = {0, 0, 0, 0};
    typedef typename accum_type<Group>::type accum;
#pragma omp parallel for
    for(lint x0 = 0; x0 < h.U->getLt(); x0++) {
      for(lint x1 = 0; x1 < h.U->getLx(); x1++) {
        for(lint x2 = 0; x2 < h.U->getLy(); x2++) {
          for(lint x3 = 0; x3 < h.U->getLz(); x3++) {
            std::vector<lint> x = {x0, x1, x2, x3};
            for(size_t mu = 0; mu < h.U->getndims(); mu++) {
              accum S;
              get_staples(S, *h.U, x, mu);
              S = (*h.U)(x, mu) * S;
              // the antihermitian traceless part
              // beta/N_c *(U*U^stap - (U*U^stap)^dagger)
              // in get_deriv
              deriv(x, mu) += fac*h.U->getBeta()/double(h.U->getNc()) * get_deriv<double>(S);
            }
          }
        }
      }
    }
    return;
  }
};
