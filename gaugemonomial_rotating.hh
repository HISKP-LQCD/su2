// gaugemonomial_rotating.hh
/**
 * @file gaugemonomial_rotating.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief
 * @version 0.1
 * @date 2022-05-05
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include "accum_type.hh"
#include "adjointfield.hh"
#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "hamiltonian_field.hh"
#include "monomial.hh"
#include "su2.hh"
#include "u1.hh"
#include <complex>
#include <vector>

namespace rotating_frame {
  template <class T> using nd_max_arr = spacetime_lattice::nd_max_arr<T>;

  /**
   * @brief return x+\hat{\mu}
   * @param x
   * @param mu
   * @return nd_max_arr
   */
  nd_max_arr<size_t> xp(nd_max_arr<size_t> x, const size_t &mu) {
    x[mu]++; // note: 'x' is passed value (this is a copy)
    return x;
  }

  /**
   * @brief return x-\hat{\mu}
   * @param x
   * @param mu
   * @return nd_max_arr
   */
  nd_max_arr<size_t> xm(nd_max_arr<size_t> x, const size_t &mu) {
    x[mu]--; // note: 'x' passed by value (this is a copy)
    return x;
  }

  /**
   * @brief Get the plaquette U_{\mu\nu} as in eq. (2.48) of
   * https://link.springer.com/book/10.1007/978-3-642-01850-3
   * @tparam Group
   * @param U gauge links
   * @param x
   * @param mu
   * @param nu
   * @return double
   */
  template <class Group>
  Group plaquette(const gaugeconfig<Group> &U,
                  const nd_max_arr<size_t> &x,
                  const size_t &mu,
                  const size_t &nu) {
    const Group res =
      U(x, mu) * U(xp(x, mu), nu) * U(xp(x, nu), mu).dagger() * U(x, nu).dagger();
    return res;
  }

  /**
   * @brief Get the clover leaf plaquette Ubar as in eq. 16 of
   * https://arxiv.org/pdf/1303.6292.pdf See also eq. (9.12) of
   * https://link.springer.com/book/10.1007/978-3-642-01850-3
   *
   * Note: In the end, we always sum over all points of the lattice.
   * Hence, in flat spacetime we can replace this clover average at each point with the
   * standard plaquette. Conversely, in curved space-time we can't because the
   * coefficients depend on the spacetime coordinates.
   * (cit. https://arxiv.org/pdf/1303.6292.pdf)
   * @tparam Group
   * @param U gauge links
   * @param x
   * @param mu
   * @param nu
   * @return double
   */
  template <class Group>
  double retr_clover_leaf(const gaugeconfig<Group> &U,
                          const nd_max_arr<size_t> &x,
                          const size_t &mu,
                          const size_t &nu) {
    double res = 0.0;

    res += retrace(plaquette(U, x, mu, nu));
    res += retrace(plaquette(U, xm(x, mu), mu, nu));
    res += retrace(plaquette(U, xm(xm(x, mu), nu), mu, nu));
    res += retrace(plaquette(U, xm(x, nu), mu, nu));

    return res / 4.0;
  }

  /**
   * @brief chair loop
   * chair loop built as the product of 2 orthogonal plaquettes: see eq. 17 of
   * https://arxiv.org/pdf/1303.6292.pdf
   * @tparam Group
   * @param U
   * @param x_munu origin of the 1st plaquette (plane mu, nu)
   * @param x_nurho origin of the 2nd plaquette (plane nu, rho)
   * @param mu
   * @param nu
   * @param rho
   * @return double
   */
  template <class Group>
  double retr_chair(const gaugeconfig<Group> &U,
                    const nd_max_arr<size_t> &x_munu,
                    const nd_max_arr<size_t> &x_nurho,
                    const size_t &mu,
                    const size_t &nu,
                    const size_t &rho) {
    return retrace(plaquette(U, x_munu, mu, nu) * plaquette(U, x_nurho, nu, rho));
  }

  /**
   * @brief Get the first chair loop: first term in the parentheses of eq. 17 of
   * https://arxiv.org/pdf/1303.6292.pdf
   * @tparam Group
   * @param U gauge config pointer
   * @param x center of the loop
   * @param mu
   * @param nu
   * @param rho
   * @return double
   */
  template <class Group>
  double retr_chair_1(const gaugeconfig<Group> &U,
                      const nd_max_arr<size_t> &x,
                      const size_t &mu,
                      const size_t &nu,
                      const size_t &rho) {
    double res = 0.0;

    res += retr_chair(U, x, x, mu, nu, rho);
    res += retr_chair(U, xm(x, nu), xm(x, nu), mu, nu, rho);

    res += retr_chair(U, xm(x, mu), xm(x, rho), mu, nu, rho);
    res += retr_chair(U, xm(xm(x, mu), nu), xm(xm(x, rho), nu), mu, nu, rho);

    return res;
  }

  /**
   * @brief Get the second chair loop: second term in the parentheses of eq. 17 of
   * https://arxiv.org/pdf/1303.6292.pdf
   * @tparam Group
   * @param U gauge config pointer
   * @param x center of the loop
   * @param mu
   * @param nu
   * @param rho
   * @return double
   */
  template <class Group>
  double retr_chair_2(const gaugeconfig<Group> &U,
                      const nd_max_arr<size_t> &x,
                      const size_t &mu,
                      const size_t &nu,
                      const size_t &rho) {
    double res = 0.0;

    res += retr_chair(U, x, xm(x, rho), mu, nu, rho);
    res += retr_chair(U, xm(x, nu), xm(xm(x, rho), nu), mu, nu, rho);

    res += retr_chair(U, xm(x, mu), x, mu, nu, rho);
    res += retr_chair(U, xm(xm(x, mu), nu), xm(x, nu), mu, nu, rho);

    return res;
  }

  /**
   * @brief get the anti-symmetric average of the chair loop: eq. (17) of eq. 17 of
   * https://arxiv.org/pdf/1303.6292.pdf
   *
   * @tparam Group
   * @param U
   * @param x
   * @param mu
   * @param nu
   * @param rho
   * @return double
   */
  template <class Group>
  double retr_asymm_chair(const gaugeconfig<Group> &U,
                          const nd_max_arr<size_t> &x,
                          const size_t &mu,
                          const size_t &nu,
                          const size_t &rho) {
    return (retr_chair_1(U, x, mu, nu, rho) - retr_chair_2(U, x, mu, nu, rho)) / 8.0;
  }

  /**
   * @brief analogous of gauge energy in flat spacetime (Minkowsky metric)
   * sum of plaquettes such that in the limit Omega=0.0 one obtains the usual gauge energy
   * (see `su2/gauge_energy.hh`)
   * @tparam Group
   * @param U
   * @param Omega
   * @return double
   */
  template <class Group>
  double gauge_energy(const gaugeconfig<Group> &U,
                      const double &Omega,
                      const bool &spatial_only = false) {
    const double Omega2 = std::pow(Omega, 2);

    // 0 if spatial_only==false, 1 if spatial_only==true
    size_t startmu = (size_t)spatial_only;
    double res = 0.0;
#pragma omp parallel for reduction(+ : res)
    for (size_t x0 = 0; x0 < U.getLt(); x0++) {
      for (size_t x1 = 0; x1 < U.getLx(); x1++) {
        for (size_t x2 = 0; x2 < U.getLy(); x2++) {
          for (size_t x3 = 0; x3 < U.getLz(); x3++) {
            const nd_max_arr<size_t> x = {x0, x1, x2, x3};
            const double r2 = x1 * x1 + x2 * x2;

            res += retr_clover_leaf(U, x, 0, 1);
            res += retr_clover_leaf(U, x, 0, 2);
            res += retr_clover_leaf(U, x, 0, 3);
            res += (1 + r2 * Omega2) * retr_clover_leaf(U, x, 1, 2);
            res += (1 + x2 * x2 * Omega2) * retr_clover_leaf(U, x, 1, 3);
            res += (1 + x1 * x1 * Omega2) * retr_clover_leaf(U, x, 2, 3);

            res += x2 * Omega * retr_asymm_chair(U, x, 0, 1, 2);
            res -= x1 * Omega * retr_asymm_chair(U, x, 0, 2, 1);
            res += x2 * Omega * retr_asymm_chair(U, x, 0, 1, 3);
            res -= x1 * Omega * retr_asymm_chair(U, x, 0, 2, 3);
            res += x1 * x2 * Omega2 * retr_asymm_chair(U, x, 1, 3, 2);
          }
        }
      }
    }
    return res;
  }

  template <class Group>
  double get_S_G(const gaugeconfig<Group> &U, const double &Omega) {
    return U.getBeta() *
           (U.getVolume() * 6 - gauge_energy<Group>(U, Omega) / double(U.getNc()));
  }

  /**
   * @brief gauge force in the HMC
   *
   * @tparam Group
   * @param U
   * @param Omega
   * @return Group
   */


template<class T, class G> T get_F_G(const gaugeconfig<G> &U,
                const nd_max_arr<size_t> &x,
                const size_t &mu,
                const double &Omega) {
    const T Stap = (*U)(x, mu) * get_staples(U, x, mu);
    return U.getBeta() / double(U.getNc()) * get_deriv<double>(Stap);
  }

  /**
   * @brief gauge monomial in a rotating frame of reference
   * class describing the gauge monomial S_G in a rotating frame of reference as in
   * https://arxiv.org/pdf/1303.6292.pdf Without loss of generality, it is always assumed
   * that the rotation is around the 'z' axis.
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
     * @brief derivative with thespect to the gauge field
     * Notes:
     *   1. In the action S_G, for each point 'x' the staples come from the plaquette at
     * the point itself, plus the contribution from nearest neighbors. All of them loop
     * "clock-wise". However, since in the end we take the Real part of the Trace, we can
     * replace the latter loops by their hermitian conjugate --> counterclockwise, and
     * also put U_{\mu}(x) in front. Therefore, in flat spacetime: S_G = (\beta)* \sum_{x}
     * [1 - (1/N_c)*Re[ w(x) Tr[ U_\mu{x}*S_{\mu}(x) ] ],
     *   2. In flat spacetime the action contains the sum of all plaquettes,
     *      so the staple contains 1 contribution from the plaquette at 'x',
     *      and (nd-1) from the others. All of them have the same "spacetime" weight (flat
     * metric)
     *   3. In curved spacetime we consider the clover leaf plaquette,
     *      hence there are 2*(nd-1) contributions with the "spacetime" weight at 'x'
     *      and other 2*(nd-1) weighted with the metric at the nearest neighbors points.
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

              // deriv(x, mu) += fac * get_F_G<accum>(*h.U, x, mu, Omega);

              accum S;
              get_staples(S, *h.U, x, mu);
              S = (*h.U)(x, mu) * S;              
              deriv(x, mu) += fac*h.U->getBeta()/double(h.U->getNc()) * get_deriv<double>(S);
              
              }
            }
          }
        }
      }
      return;
    }
  };

} // namespace rotating_frame
