/**
 * @file gaugemonomial_rotating.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief
 * @version 0.1
 * @date 2022-05-05
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include "rotating-gauge_energy.hpp"

namespace rotating_spacetime {

  namespace rot_st = rotating_spacetime;

  /**
   * @brief gauge action S_G as in eq. 15 of https://arxiv.org/pdf/1303.6292.pdf
   * @tparam Group gauge symmetry group
   * @param U gauge configuration
   * @param Omega imaginary angular frequency
   * @return double value of S_G
   */
  template <class Group>
  double get_S_G(const gaugeconfig<Group> &U, const double &Omega) {
    const double Lt = U.getLt(), Lx = U.getLx(), Ly = U.getLy(), Lz = U.getLz();
    const double ndims_fact = spacetime_lattice::num_pLloops_half(U.getndims());

    const double sum_x2 = Lx * (Lx + 1) * (2 * Lx + 1) / 6.0; // \sum_{x_1=1}^{Lx} x_1^2
    const double sum_y2 = Ly * (Ly + 1) * (2 * Ly + 1) / 6.0; // \sum_{x_2=1}^{Ly} x_2^2

    const double one = U.getVolume() * ndims_fact + Lt * Lz * (U.getndims() - 2) *
                                                      (sum_x2 + sum_y2) *
                                                      std::pow(Omega, 2);
    const double two = rot_st::gauge_energy<Group>(U, Omega) / double(U.getNc());

    return U.getBeta() * (one - two);
  }


  /**
   * @brief gauge force in the HMC
   *
   * @tparam Group
   * @param U
   * @param Omega
   * @return Group
   */
  template <class T, class G>
  T get_F_G(const gaugeconfig<G> &U,
            const nd_max_arr<size_t> &x,
            const size_t &mu,
            const double &Omega) {
    T St1;

    for (size_t nu = 0; nu < U.getndims(); nu++) {
      if (mu == nu) {
        continue;
      }

      T S_tt = get_staple<T, G>(U, x, mu, nu, true, true);
      T S_tf = get_staple<T, G>(U, x, mu, nu, true, false);
      T S_ft = get_staple<T, G>(U, x, mu, nu, false, true);
      T S_ff = get_staple<T, G>(U, x, mu, nu, false, false);

      // the 2 contributions from the clover plaquette at x
      St1 += (S_tt + S_ff); // plaq_factor(x, mu, nu, Omega) * (S_tt + S_ff);

      /* contributions from nearest neighbors */

      // x+\nu
      St1 += plaq_factor(xp(x, nu), mu, nu, Omega) * S_tt;
      // x+\mu+\nu
      St1 += plaq_factor(xp(xp(x, mu), nu), mu, nu, Omega) * S_tt;
      // x+\mu
      St1 += plaq_factor(xp(x, mu), mu, nu, Omega) * (S_tt + S_ff);
      // x-\nu
      St1 += plaq_factor(xm(x, nu), mu, nu, Omega) * S_ff;
      // x+\mu-\nu
      St1 += plaq_factor(xp(xm(x, nu), mu), mu, nu, Omega) * S_ff;
    }
    St1 /= 4.0; // each plaquette comes from an average over the clover

    T St2;

    // Stap += x2 * Omega * retr_asymm_chair(U, x, 0, 1, 2);
    // Stap -= x1 * Omega * retr_asymm_chair(U, x, 0, 2, 1);
    // Stap += x2 * Omega * retr_asymm_chair(U, x, 0, 1, 3);
    // Stap -= x1 * Omega * retr_asymm_chair(U, x, 0, 2, 3);
    // Stap += x1 * x2 * Omega2 * retr_asymm_chair(U, x, 1, 3, 2);

    T Stap = U(x, mu) * (St1 + St2);

    return Stap;
  }

  /**
   * @brief gauge monomial in a rotating frame of reference
   * class describing the gauge monomial S_G in a rotating frame of reference as in
   * https://arxiv.org/pdf/1303.6292.pdf.
   * Without loss of generality, it is assumed a rotation around the 'z' axis.
   * @tparam Float
   * @tparam Group
   */
  template <typename Float, class Group>
  class gauge_monomial : public monomial<Float, Group> {
  private:
    double Omega; // imaginary angular velocity
  public:
    gauge_monomial<Float, Group>(unsigned int _timescale, const double &_Omega)
      : monomial<Float, Group>::monomial(_timescale), Omega(_Omega) {}

    void heatbath(hamiltonian_field<Float, Group> const &h) override {
      monomial<Float, Group>::Hold = get_S_G(*h.U, (*this).Omega);
      return;
    }

    void accept(hamiltonian_field<Float, Group> const &h) override {
      monomial<Float, Group>::Hnew = get_S_G(*h.U, (*this).Omega);
      return;
    }

    /**
     * @brief derivative with respect to the gauge field
     * Notes:
     *   1. In the action S_G, for each point 'x' the staples come from the plaquette at
     * the point itself, plus the contribution from nearest neighbors. All of them loop
     * "counterclockwise". However, since in the end we take the Real part of the Trace,
     * we can replace the latter loops by their hermitian conjugate --> clockwise, and
     * also put U_{\mu}(x) in front. Therefore, in flat spacetime: S_G = (\beta)* \sum_{x}
     * [1 - (1/N_c)*Re[ w(x) Tr[ U_\mu{x}*S_{\mu}(x) ] ],
     *   2. In flat spacetime the action contains the sum of all plaquettes,
     *      so the staple contains 1 contribution from the plaquette at 'x',
     *      and (nd-1) from the others. All of them have the same "spacetime" weight (flat
     * metric)
     *   3. In curved spacetime we consider the clover leaf plaquette,
     *      hence there are 2*(nd-1) contributions with the "spacetime" weight at 'x'
     *      and other 6*(nd-1) weighted with the metric at the nearest neighbors points.
     * @param deriv reference to the derivative object
     * @param h hamiltonian field
     * @param fac
     */
    void derivative(adjointfield<Float, Group> &deriv,
                    hamiltonian_field<Float, Group> const &h,
                    const Float fac = 1.) const override {
      typedef typename accum_type<Group>::type accum;
#pragma omp parallel for
      for (size_t x0 = 0; x0 < h.U->getLt(); x0++) {
        for (size_t x1 = 0; x1 < h.U->getLx(); x1++) {
          for (size_t x2 = 0; x2 < h.U->getLy(); x2++) {
            for (size_t x3 = 0; x3 < h.U->getLz(); x3++) {
              const nd_max_arr<size_t> x = {x0, x1, x2, x3};
              for (size_t mu = 0; mu < h.U->getndims(); mu++) {
                accum F = get_F_G<accum>(*h.U, x, mu, Omega);
                deriv(x, mu) +=
                  fac * h.U->getBeta() / double(h.U->getNc()) * get_deriv<double>(F);
              }
            }
          }
        }
      }
      return;
    }
  };

} // namespace rotating_spacetime
