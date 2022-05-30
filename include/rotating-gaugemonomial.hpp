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
    const double Omega; // imaginary angular velocity
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
                accum F = get_staples_with_st_fact<accum>(*h.U, x, mu, Omega);
                deriv(x, mu) += fac * get_deriv<double>(F);
              }
            }
          }
        }
      }
      return;
    }
  };

} // namespace rotating_spacetime
