/**
 * @file flat_spacetime_gauge_energy.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief gauge energy in flat spacetime (euclidean metric)
 * @version 0.1
 * @date 2022-05-20
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

#ifdef _USE_OMP_
#include <omp.h>
#endif

namespace rotating_spacetime {
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
   * @param U gauge configuration
   * @param x_munu origin of the 1st plaquette (plane mu, nu)
   * @param x_nurho origin of the 2nd plaquette (plane nu, rho)
   * @param mu
   * @param nu
   * @param rho
   * @param flip  true if orientation is flipped for the 2nd plaquette
   * @return double
   */
  template <class Group>
  Group chair_loop(const gaugeconfig<Group> &U,
                   const nd_max_arr<size_t> &x_munu,
                   const nd_max_arr<size_t> &x_nurho,
                   const size_t &mu,
                   const size_t &nu,
                   const size_t &rho,
                   const bool &flip = true) {
    Group C = plaquette(U, x_munu, mu, nu);

    if (flip) {
      C *= plaquette(U, x_nurho, nu, rho).dagger();
    } else {
      C *= plaquette(U, x_nurho, nu, rho);
    }

    return C;
  }

  // Re(Tr(chair))
  template <class Group>
  double retr_chair_loop(const gaugeconfig<Group> &U,
                         const nd_max_arr<size_t> &x_munu,
                         const nd_max_arr<size_t> &x_nurho,
                         const size_t &mu,
                         const size_t &nu,
                         const size_t &rho,
                         const bool &flip = true) {
    return retrace(chair_loop(U, x_munu, x_nurho, mu, nu, rho, flip));
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
  double retr_chair_loop_1(const gaugeconfig<Group> &U,
                           const nd_max_arr<size_t> &x,
                           const size_t &mu,
                           const size_t &nu,
                           const size_t &rho) {
    double res = 0.0;

    res += retr_chair_loop(U, x, x, mu, nu, rho);
    res += retr_chair_loop(U, xm(x, nu), xm(x, nu), mu, nu, rho);

    res += retr_chair_loop(U, xm(x, mu), xm(x, rho), mu, nu, rho);
    res += retr_chair_loop(U, xm(xm(x, mu), nu), xm(xm(x, rho), nu), mu, nu, rho);

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
  double retr_chair_loop_2(const gaugeconfig<Group> &U,
                           const nd_max_arr<size_t> &x,
                           const size_t &mu,
                           const size_t &nu,
                           const size_t &rho) {
    double res = 0.0;

    res += retr_chair_loop(U, x, xm(x, rho), mu, nu, rho, true);
    res += retr_chair_loop(U, xm(x, nu), xm(xm(x, rho), nu), mu, nu, rho, true);

    res += retr_chair_loop(U, xm(x, mu), x, mu, nu, rho, true);
    res += retr_chair_loop(U, xm(xm(x, mu), nu), xm(x, nu), mu, nu, rho, true);

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
    return (retr_chair_loop_1(U, x, mu, nu, rho) - retr_chair_loop_2(U, x, mu, nu, rho)) /
           8.0;
  }

  /**
   * @brief analogous of gauge energy in flat spacetime (Minkowsky metric)
   * sum of plaquettes such that in the limit Omega=0.0 one obtains the usual gauge energy
   * (see `gauge_energy.hpp`)
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

} // namespace rotating_spacetime