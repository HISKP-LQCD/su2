#pragma once
//#include"gradient_flow.hh"
#include "adjointfield.hh"
#include "flat-energy_density.hh"
#include "flat-gauge_energy.hpp"
#include "gaugeconfig.hh"
#include "flat-gaugemonomial.hh"
#include "hamiltonian_field.hh"
#include "monomial.hh"
#include "su2.hh"
#include "update_gauge.hh"

#include <fstream>
#include <iomanip>
#include <string>

template <typename Float, class Group>
void runge_kutta(hamiltonian_field<Float, Group> &h,
                 monomial<Float, Group> &SW,
                 const double eps) {
  double zfac[5] = {(-17.0) / (36.0), (8.0) / (9.0), (-3.0) / (4.0)};
  double expfac[3] = {-36.0 / 4. / 17.0, 1., -1.};

  // following arxiv:1006.4518
  //
  // W0      = V_t
  // W1      = exp(1/4 Z0) W0
  // W2      = exp(8/9 Z1 - 17/36 Z0) W1 = exp(Z') W1
  // V_t+eps = exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) W2 = exp(3/4 Z2 - Z') W2
  //
  // Zi = eps*Z(Wi)
  // before the three steps zero the derivative field
  zeroadjointfield(*(h.momenta));

  for (int f = 0; f < 3; f++) {
    // add to *(h.momenta)
    // we have to cancel beta/N_c from the derivative
    // a factor two to obtain the correct normalisation of
    // the Wilson plaquette action
    // S_W = 1./g_0^2 \sum_x \sum_{p} Re Tr(1 - U(p))
    // where we sum over all oriented plaquettes
    // we sum over unoriented plaquettes, so we have to multiply by 2
    // which is usually in beta
    SW.derivative(*(h.momenta), h, 2. * h.U->getNc() * zfac[f] / h.U->getBeta());
    // The '-' comes from the action to be tr(1-U(p))
    // update the flowed gauge field Vt
    update_gauge(h, -eps * expfac[f]);
  }
  return;
}

namespace flat_spacetime {

  template <class Group>
  void gradient_flow(const gaugeconfig<Group> &U,
                     std::string const &path,
                     const double tmax,
                     const double &eps) {
    const double ndims_fact = spacetime_lattice::num_pLloops_half(U.getndims());

    double t[3];
    double P[3], E[3], Q[3];
    double P_ss[3], E_ss[3], Q_ss[3]; // spatial-spatial contribution

    std::ofstream os(path, std::ios::out);

    double density = 0., topQ = 0.; // dummy variables
    double density_ss = 0., topQ_ss = 0.; // dummy variables
    for (unsigned int i = 0; i < 3; i++) {
      t[i] = 0.;
      P[i] = 0.;
      E[i] = 0.;
      Q[i] = 0.;
      P_ss[i] = 0.;
      E_ss[i] = 0.;
      Q_ss[i] = 0.;
    }
    const double den = U.getVolume() * double(U.getNc()) * ndims_fact;
    P[2] = flat_spacetime::gauge_energy(U) / den;
    P_ss[2] = flat_spacetime::gauge_energy(U, true) / den;
    flat_spacetime::energy_density(U, density, topQ);
    E[2] = density;
    flat_spacetime::energy_density(U, density_ss, topQ_ss, true, true);
    E_ss[2] = density_ss;

    // definine a fictitious gauge configuration Vt, momenta and hamiltonian field to
    // evolve with the flow
    gaugeconfig<Group> Vt(U);
    adjointfield<double, Group> deriv(U.getLx(), U.getLy(), U.getLz(), U.getLt(),
                                      U.getndims());
    hamiltonian_field<double, Group> h(deriv, Vt);
    gaugemonomial<double, Group> SW(
      0); // gradient flow for the Wilson (pure) gauge action

    os << "t "; // flow time
    os << "P P_ss "; // plaquette
    os << "Ep Ep_ss "; // energy from regular plaquette
    os << "t^2*Ep(t) t^2*Ep_ss(t) "; // t^2 * E_p(t)
    os << "Ec Ec_ss "; // energy from clover-leaf plaquette
    os << "t^2*Ec(t) t^2*Ec_ss(t) "; // t^2 * E_c(t)
    os << "Q Q_ss" << std::endl;

    // evolution of t[1] until tmax
    //(note: eps=0.01 and tmax>0 --> the loop ends at some point)
    // at each step we consider a triplet of values for t,P,E,Q
    while (t[1] < tmax) {
      // splicing the results at t[2] to the new o-th temporal time slice
      t[0] = t[2];
      P[0] = P[2];
      E[0] = E[2];
      Q[0] = Q[2];
      P_ss[0] = P_ss[2];
      E_ss[0] = E_ss[2];
      Q_ss[0] = Q_ss[2];

      for (unsigned int x0 = 1; x0 < 3; x0++) {
        t[x0] = t[x0 - 1] + eps;
        runge_kutta(h, SW, eps); //
        P[x0] = flat_spacetime::gauge_energy(Vt) / den;
        P_ss[x0] = flat_spacetime::gauge_energy(Vt, true) / den;
        flat_spacetime::energy_density(Vt, density, topQ);
        E[x0] = density;
        Q[x0] = topQ;
        flat_spacetime::energy_density(Vt, density_ss, topQ_ss, true);
        E_ss[x0] = density_ss;
        Q_ss[x0] = topQ_ss;
      }

      const double tsqr = t[1] * t[1];

      const double factP = 2 * U.getNc() * ndims_fact; // normalization factor

      const double Ep = factP * (1. - P[1]);
      const double Ep_ss = factP * (1. - P_ss[1]);

      const double t2Ep = tsqr * Ep;
      const double t2Ep_ss = tsqr * Ep_ss;

      const double Ec = E[1];
      const double Ec_ss = E_ss[1];

      const double t2Ec = tsqr * Ec;
      const double t2Ec_ss = tsqr * Ec_ss;

      os << std::scientific; // using scientifc notation
      os.precision(16);
      os << t[1] << " ";
      os << P[1] << " " << P_ss[1] << " ";
      os << Ep << " " << Ep_ss << " ";
      os << t2Ep << " " << t2Ep_ss << " ";
      os << Ec << " " << Ec_ss << " ";
      os << t2Ec << " " << t2Ec_ss << " ";
      os << Q[1] << " " << Q_ss[1] << "\n";
    }

    return;
  }

} // namespace flat_spacetime
