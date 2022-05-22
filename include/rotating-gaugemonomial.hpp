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

  template <class Group>
  double get_S_G(const gaugeconfig<Group> &U, const double &Omega) {
    return U.getBeta() *
           (U.getVolume() * 6 -
            rotating_spacetime::gauge_energy<Group>(U, Omega) / double(U.getNc()));
  }

  /**
   * @brief Get the staple object S_{\mu \nu}(x)
   * return S_{\mu \nu}(x) = U_{\nu}(x+\mu) U_{\mu}^{\dagger}(x+\nu) U_{\nu}^{\dagger}(x)
   * in the adjoint representation T
   * @tparam T --> has to be the adjoint representation
   * @tparam S gauge group
   * @tparam Arr
   * @param U gauge configuration
   * @param x point of application of the staple
   * @param mu
   * @param nu
   * @param clockwise  true when the staple orientation is clock-wise
   */
  template <class T, class S, class Arr>
  T get_staple(const gaugeconfig<S> &U,
               const Arr &x,
               const size_t &mu,
               const size_t &nu,
               const bool &clockwise) {
    T K;
    if (clockwise) {
      K = U(xp(x, mu), nu) * U(xp(x, nu), mu).dagger() * U(x, nu).dagger();
    } else {
      K = U(xm(xp(x, mu), nu), nu).dagger() * U(xm(x, nu), mu).dagger() * U(x, nu);
    }
    return K;
  }

  /**
   * @brief spacetime metric factor multiplying the staple
   * @param mu
   * @param nu
   * @return double
   */
  double staple_factor(const nd_max_arr<size_t> &x,
                       const size_t &mu,
                       const size_t &nu,
                       const double &Omega) {
    if (mu == 0) {
      return 1.0;
    }
    if (mu == 1 && nu == 2) {
      return (1 + (std::pow(x[1], 2) + std::pow(x[2], 2)) * std::pow(Omega, 2));
    }
    if (mu == 1 && nu == 3) {
      return (1 + std::pow(x[2], 2) * std::pow(Omega, 2));
    }
    if (mu == 2 && nu == 3) {
      return (1 + std::pow(x[1], 2) * std::pow(Omega, 2));
    } else {
#pragma omp single
      {
        std::cout << "Error: mu=" << mu << " and nu=" << nu << " not supported by "
                  << __func__ << "\n";
      }
      abort();
      return 0.0;
    }
  }

  /**
   * @brief Get the chair staple object S_{\mu \nu \rho}(x)
   * return S_{\mu \nu \rho}(x) = U_{\mu}(x)^{\dagger}*chair_loop,
   * where chair_loop is product of 2 orthogonal plaquettes: see eq. 17 of
   * https://arxiv.org/pdf/1303.6292.pdf.
   * The result is computed as the product of a staple and a plaquette
   * @tparam T --> has to be the adjoint representation
   * @tparam S gauge group
   * @tparam Arr
   * @param U gauge configuration
   * @param x point of application of the staple
   * @param mu
   * @param nu
   * @param clockwise1  true when the staple orientation is clock-wise
   * @param clockwise2  true when plaquette orientation is clock-wise
   */
  template <class T, class S, class Arr>
  T get_chair_staple(gaugeconfig<S> &U,
                     const Arr &x,
                     const size_t &mu,
                     const size_t &nu,
                     const size_t &rho,
                     const bool &clockwise1,
                     const bool &clockwise2) {
    T K;
    if (clockwise1 && clockwise2) {
      K = get_staple(U, x, mu, nu, clockwise1) *
          rotating_spacetime::plaquette(U, x, nu, rho);
    }
    if (clockwise1 && !clockwise2) {
      K = rotating_spacetime::plaquette(U, xm(x, nu), nu, rho);
    }
    if (!clockwise1 && clockwise2) {
      K = rotating_spacetime::plaquette(U, xm(x, rho), nu, rho);
    }
    if (!clockwise1 && !clockwise2) {
      K = rotating_spacetime::plaquette(U, xm(xm(x, nu), rho), nu, rho);
    }
    return K;
  }

  /**
   * @brief spacetime metric factor multiplying the chair staple
   * @param mu
   * @param nu
   * @return double
   */
  double chair_staple_factor(const nd_max_arr<size_t> &x,
                             const size_t &mu,
                             const size_t &nu,
                             const size_t &rho,
                             const double &Omega) {
    if (mu == 0) {
      if (nu == 1 && rho == 2) {
        return x[2] * Omega;
      }
      if (nu == 2 && rho == 1) {
        return -x[1] * Omega;
      }
      if (nu == 1 && rho == 3) {
        return x[2] * Omega;
      }
      if (nu == 2 && rho == 3) {
        return -x[1] * Omega;
      }
    }
    if (mu == 1 && nu == 3 && rho == 2) {
      return x[1] * x[2] * std::pow(Omega, 2);
    } else {
      std::cout << "Error: mu=" << mu << " and nu=" << nu << " not supported by "
                << __func__ << "\n";
      abort();
      return 0.0;
    }
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

    for (size_t nu = mu + 1; nu < U.getndims(); nu++) {
      // the 2 contributions from the clover plaquette at x
      St1 += staple_factor(x, mu, nu, Omega) * (get_staple<T, G>(U, x, mu, nu, true) +
                                                get_staple<T, G>(U, x, mu, nu, false));

      /* contributions from nearest neighbors */

      // x+nu
      St1 +=
        staple_factor(xp(x, nu), mu, nu, Omega) * get_staple<T, G>(U, x, mu, nu, true);
      // at x-nu
      St1 +=
        staple_factor(xm(x, nu), mu, nu, Omega) * get_staple<T, G>(U, x, mu, nu, false);
      //  x+mu+nu
      St1 += staple_factor(xp(xp(x, mu), nu), mu, nu, Omega) *
             get_staple<T, G>(U, x, mu, nu, true);
      //  x-nu
      St1 += staple_factor(xp(xm(x, nu), mu), mu, nu, Omega) *
             get_staple<T, G>(U, x, mu, nu, false);
    }

    St1 /= 4.0; // each plaquette comes from an average over the clover

    T St2;

    // Stap += x2 * Omega * retr_asymm_chair(U, x, 0, 1, 2);
    // Stap -= x1 * Omega * retr_asymm_chair(U, x, 0, 2, 1);
    // Stap += x2 * Omega * retr_asymm_chair(U, x, 0, 1, 3);
    // Stap -= x1 * Omega * retr_asymm_chair(U, x, 0, 2, 3);
    // Stap += x1 * x2 * Omega2 * retr_asymm_chair(U, x, 1, 3, 2);

    T Stap = U(x, mu) * (St1 + St2);
    return U.getBeta() / double(U.getNc()) * Stap;
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
                accum F = get_F_G<accum>(*h.U, x, mu, Omega);
                deriv(x, mu) +=
                  fac * h.U->getBeta() / double(h.U->getNc()) * get_deriv<double>(F);

                // accum S;
                // get_staples(S, *h.U, x, mu);
                // S = (*h.U)(x, mu) * S;
                // deriv(x, mu) += fac*h.U->getBeta()/double(h.U->getNc()) *
                // get_deriv<double>(S);
              }
            }
          }
        }
      }
      return;
    }
  };

} // namespace rotating_spacetime
