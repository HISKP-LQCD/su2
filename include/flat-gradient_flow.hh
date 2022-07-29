#pragma once
//#include"gradient_flow.hh"
#include "adjointfield.hh"
#include "flat-energy_density.hh"
#include "flat-gauge_energy.hpp"
#include "flat-gaugemonomial.hh"
#include "gaugeconfig.hh"
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

  /**
   * @brief Printing Wilson gradient flow evolution on 'path'
   * Prints a dataframe for the Wilson flow evolution of
   * flowtime, Plaquette(full, spatial-spatial, temporal-spatial), Energy_plaquette(...),
   * Energy_clover[improved formula](...), etc.
   * Notes:
   * - Each type of plaquette has limit 1 --> different normalization factors for P_ss and
   * P_ts.
   * - E_i = 1 - P_i, i="","ss","ts"
   * - E = E_ss + E_ts
   * @tparam Group
   * @param U
   * @param path
   * @param tmax
   * @param eps
   */
  template <class Group>
  void gradient_flow(const gaugeconfig<Group> &U,
                     std::string const &path,
                     const double &tmax,
                     const double &eps,
                     const double &xi = 1.0) {
    const size_t d = U.getndims();
    const double ndims_fact = spacetime_lattice::num_pLloops_half(d);
    const double ndims_fact_ss = spacetime_lattice::num_pLloops_half(d - 1);

    double t[3];
    double P[3], E[3], Q[3];
    double P_ss[3], E_ss[3], Q_ss[3]; // spatial-spatial contribution

    std::ostringstream oss;

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

    const double den_common = U.getVolume() * double(U.getNc());
    // const double den = den_common * ndims_fact;
    // const double den_ss = den_common * ndims_fact_ss;
    P[2] = flat_spacetime::gauge_energy(U) / den_common;
    P_ss[2] = flat_spacetime::gauge_energy(U, true) / den_common;
    flat_spacetime::energy_density(U, density, topQ);
    E[2] = density;
    flat_spacetime::energy_density(U, density_ss, topQ_ss, true, true);
    E_ss[2] = density_ss;

    // definine a fictitious gauge configuration Vt, momenta and hamiltonian field to
    // evolve with the flow
    gaugeconfig<Group> Vt(U);
    adjointfield<double, Group> deriv(U.getLx(), U.getLy(), U.getLz(), U.getLt(), d);
    hamiltonian_field<double, Group> h(deriv, Vt);
    gaugemonomial<double, Group> SW(0, xi); // Wilson (pure) gauge action

    oss << "t "; // flow time
    oss << "P P_ss P_ts "; // plaquette
    oss << "Ep Ep_ss Ep_ts "; // energy from regular plaquette
    oss << "xi "; // t^2 * E_plaquette(t)
    oss << "Ec Ec_ss Ec_ts "; // energy from clover-leaf plaquette
    oss << "Q Q_ss Q_ts" << std::endl;

    // evolution of t[1] until tmax
    //(note: eps=0.01 and tmax>0 --> the loop ends at some point)
    // at each step we consider a triplet of values for t,P,E,Q

    while (t[1] < tmax) {
      // splicing the results at t[2] to the new 0-th temporal time slice
      t[0] = t[2];
      P[0] = P[2];
      E[0] = E[2];
      Q[0] = Q[2];
      P_ss[0] = P_ss[2];
      E_ss[0] = E_ss[2];
      Q_ss[0] = Q_ss[2];

      for (unsigned int x0 = 1; x0 < 3; x0++) {
        t[x0] = t[x0 - 1] + eps;
        runge_kutta(h, SW, eps); // apply Runge-Kutta integration method
        P[x0] = flat_spacetime::gauge_energy(Vt) / den_common;
        P_ss[x0] = flat_spacetime::gauge_energy(Vt, true) / den_common;
        flat_spacetime::energy_density(Vt, density, topQ, true);
        E[x0] = density;
        Q[x0] = topQ;
        flat_spacetime::energy_density(Vt, density_ss, topQ_ss, true, true);
        E_ss[x0] = density_ss;
        Q_ss[x0] = topQ_ss;
      }

      const double P_ts = P[1] - P_ss[1];
      const double Q_ts = Q[1] - Q_ss[1];

      const double tsqr = t[1] * t[1];

      const double Ep = ndims_fact - P[1];
      const double Ep_ss = ndims_fact_ss - P_ss[1];
      const double Ep_ts = double(d - 1.0) - P_ts;
      const double xi_R = sqrt((Ep_ts / Ep_ss) / (U.getndims() - 2)); // renormalized anisotropy

      const double t2Ep = tsqr * Ep;
      const double t2Ep_ss = tsqr * Ep_ss;
      const double t2Ep_ts = tsqr * Ep_ts;

      const double Ec = E[1];
      const double Ec_ss = E_ss[1];
      const double Ec_ts = Ec - Ec_ss;

      const double t2Ec = tsqr * Ec;
      const double t2Ec_ss = tsqr * Ec_ss;
      const double t2Ec_ts = tsqr * Ec_ts;

      oss << std::scientific; // using scientific notation
      oss.precision(16);
      oss << t[1] << " ";
      oss << P[1] << " " << P_ss[1] << " " << P_ts << " ";
      oss << Ep << " " << Ep_ss << " " << Ep_ts << " ";
      oss << xi_R << " ";
      oss << Ec << " " << Ec_ss << " " << Ec_ts << " ";
      oss << Q[1] << " " << Q_ss[1] << " " << Q_ts << "\n";
    }

    std::ofstream ofs(path, std::ios::out);
    ofs << oss.str();

    return;
  }

} // namespace flat_spacetime
