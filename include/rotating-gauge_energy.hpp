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
//#include "get_staples.hh"
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
   * @brief Get the staple object S_{\mu \nu}(x)
   * returns, in the representation T the staple
   * - in the positive/negative \nu semi-plane : up==true/false
   * - counterclockwise/clockwise oriented in the : ccwise==true/false
   * @tparam T gauge group (representation T)
   * @tparam S gauge group (representation S)
   * @tparam Arr
   * @param U gauge configuration
   * @param x point of application of the staple
   * @param mu
   * @param nu
   * @param up true when the staple is in the positive '\nu' semi-plane
   * @param ccwise  true when the staple orientation is counterclockwise
   */
  template <class T, class S, class Arr>
  T get_staple(const gaugeconfig<S> &U,
               const Arr &x,
               const size_t &mu,
               const size_t &nu,
               const bool &up,
               const bool &ccwise) {
    T K;
    S res;
    if (up) {
      res = U(xp(x, mu), nu) * U(xp(x, nu), mu).dagger() * U(x, nu).dagger();
    } else {
      res = U(xm(x, nu), nu).dagger() * U(xm(x, nu), mu) * U(xm(xp(x, mu), nu), nu);
    }

    if (ccwise) {
      K = res;
    } else {
      K = res.dagger();
    }

    return K;
  }

  /**
   * @brief Get the chair staple object CS_{\mu \nu}(x)
   * returns, in the representation T the staple
   * - in the positive/negative \nu semi-plane : up==true/false
   * - counterclockwise/clockwise oriented in the : ccwise==true/false
   * @tparam T --> has to be the adjoint representation
   * @tparam S gauge group
   * @tparam Arr
   * @param U gauge configuration
   * @param x point of application of the staple
   * @param mu
   * @param nu
   * @param p array of 2 values:
   * p[0]==true : staple is in the positive \mu semi-plane
   * p[1]==true : staple is in the positive \rho semi-plane
   * @param cc if true staple is oriented counterclockwise in the plane (\nu,\rho)
   * Note that this orientation is sufficient to characterize the chair staple,
   * since the orientation of the 1st staple coincides with p[0]
   */
  template <class T, class S, class Arr>
  T get_chair_staple(const gaugeconfig<S> &U,
                     const Arr &x,
                     const size_t &mu,
                     const size_t &nu,
                     const size_t rho,
                     const std::array<bool, 2> &p,
                     const bool &cc) {
    const bool dag_all = (p[0] ^ (p[1] && cc)); // wether to dagger in the end

    // Note: if taking the dagger in the end, the orientations are reversed
    const std::array<bool, 2> dd = {(bool)(p[0] ^ dag_all), (bool)(cc ^ dag_all)};
    S K1 = get_staple<S, S>(U, x, mu, nu, p[0], dd[0]);

    Arr y = x;
    if (!p[0]) {
      y[nu]--;
    }
    if (p[0] ^ (!(p[1] ^ cc))) {
      y[mu]++;
    }
    S M = U(y, nu); // link in the middle between the staples K1 and K2

    const bool dag_U = cc;
    if (dag_U) {
      M = M.dagger();
    }

    S K2 = get_staple<S, S>(U, y, nu, rho, p[1], cc);

    S K = K1 * M * K2;
    if (dag_all) {
      K = K.dagger();
    }

    return (T)K;
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
                  const size_t &nu,
                  const bool &up,
                  const bool &ccwise) {
    Group S = get_staple<Group, Group>(U, x, mu, nu, up, ccwise);

    if (up ^ ccwise) {
      return S * U(x, mu).dagger();
    } else {
      return U(x, mu) * S;
    }
  }

  /**
   * @brief real part and trace of the plaquette
   *
   * @tparam Group
   * @param U
   * @param x
   * @param mu
   * @param nu
   * @param up
   * @param ccwise
   * @return Group
   */
  template <class Group>
  double retr_plaquette(const gaugeconfig<Group> &U,
                        const nd_max_arr<size_t> &x,
                        const size_t &mu,
                        const size_t &nu) {
    const bool or1 = true; // orientation doesn't matter when taking the Re(Tr(...))
    return retrace(plaquette(U, x, mu, nu, true, or1));
  }

  /**
   * @brief Get the clover leaf plaquette \bar{U} as in eq. 16 of
   * https://arxiv.org/pdf/1303.6292.pdf. See also eq. (9.12) of
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

    res += retr_plaquette(U, x, mu, nu);
    res += retr_plaquette(U, xm(x, mu), mu, nu);
    res += retr_plaquette(U, xm(xm(x, mu), nu), mu, nu);
    res += retr_plaquette(U, xm(x, nu), mu, nu);

    return res / 4.0;
  }

  /**
   * @brief Get the chair loop object CL_{\mu \nu \rho}(x)
   * CL_{\mu \nu \rho}(x) is product of 2 staples.
   * see eq. 17 of https://arxiv.org/pdf/1303.6292.pdf.
   * @tparam S gauge group
   * @tparam Arr
   * @param U gauge configuration
   * @param x point of application of the staple
   * @param mu
   * @param nu
   * @param p array of bool flags for the directions \mu and \rho
   * for i=\mu,\rho : p[i]==true when the chair is in the positive semi-plane of the
   * direction 'i'.
   * The chairs in the negative '\nu' semi-plane are obtained translating the chair ar
   * x-\nu
   */
  template <class S, class Arr>
  S chair_loop(const gaugeconfig<S> &U,
               const Arr &x,
               const size_t &mu,
               const size_t &nu,
               const size_t &rho,
               const std::array<bool, 2> &p,
               const bool &cc) {
    S K = get_chair_staple<S, S>(U, x, mu, nu, rho, p, cc);
    return U(x, mu) * K;
  }

  // Re(Tr(chair))
  template <class Group>
  double retr_chair_loop(const gaugeconfig<Group> &U,
                         const nd_max_arr<size_t> &x,
                         const size_t &mu,
                         const size_t &nu,
                         const size_t &rho,
                         const std::array<bool, 2> &p) {
    const bool or1 = true; // orientation doesn't matter when taking the Re(Tr(...))
    return retrace(chair_loop(U, x, mu, nu, rho, p, or1));
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

    const nd_max_arr<size_t> xmnu = xm(x, nu);

    res += retr_chair_loop(U, x, mu, nu, rho, {true, true});
    res += retr_chair_loop(U, xmnu, mu, nu, rho, {true, true});

    res += retr_chair_loop(U, x, mu, nu, rho, {false, false});
    res += retr_chair_loop(U, xmnu, mu, nu, rho, {false, false});

    return res / 8.0;
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

    const nd_max_arr<size_t> xmnu = xm(x, nu);

    res += retr_chair_loop(U, x, mu, nu, rho, {true, false});
    res += retr_chair_loop(U, xmnu, mu, nu, rho, {true, false});

    res += retr_chair_loop(U, x, mu, nu, rho, {false, true});
    res += retr_chair_loop(U, xmnu, mu, nu, rho, {false, true});

    return res / 8.0;
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
    return (retr_chair_loop_1(U, x, mu, nu, rho) - retr_chair_loop_2(U, x, mu, nu, rho));
  }

  /**
   * @brief spacetime metric factor multiplying the chair staple
   * @param mu
   * @param nu
   * @return double
   */
  double chair_stm_factor(const nd_max_arr<size_t> &x,
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
#pragma omp single
      {
        std::cout << "Error: (mu,nu,rho)=(" << mu << "," << nu << "," << rho
                  << ") not supported by " << __func__ << "\n";
      }
      abort();
      return 0.0;
    }
  }

  /**
   * @brief spacetime metric factor multiplying the plaquette U_{\mu \nu} in the action
   * @param mu
   * @param nu
   * @return double
   */
  double plaq_factor(const nd_max_arr<size_t> &x,
                     const size_t &mu,
                     const size_t &nu,
                     const double &Omega) {
    if (mu == 0 || nu == 0) {
      return 1.0;
    }
    if ((mu == 1 && nu == 2) || (mu == 2 && nu == 1)) {
      return (1.0 + (std::pow(x[1], 2) + std::pow(x[2], 2)) * std::pow(Omega, 2));
    }
    if ((mu == 1 && nu == 3) || (mu == 3 && nu == 1)) {
      return (1.0 + std::pow(x[2], 2) * std::pow(Omega, 2));
    }
    if ((mu == 2 && nu == 3) || (mu == 3 && nu == 2)) {
      return (1.0 + std::pow(x[1], 2) * std::pow(Omega, 2));
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

            for (size_t nu = 0; nu < U.getndims(); nu++) {
              for (size_t mu = 0; mu < nu; mu++) {
                res += plaq_factor(x, mu, nu, Omega) * retr_clover_leaf(U, x, mu, nu);
              }
            }

            for (size_t mu = 0; mu < U.getndims() - 2; mu++) {
              for (size_t nu = mu + 1; nu < U.getndims() - 1; nu++) {
                for (size_t rho = 1; rho < U.getndims(); rho++) {
                  if (rho == mu || rho == nu || rho == 3) {
                    continue;
                  }
                  res += chair_stm_factor(x, mu, nu, rho, Omega) *
                         retr_asymm_chair(U, x, mu, nu, rho);
                }
              }
            }
          }
        }
      }
    }
    return res;
  }

} // namespace rotating_spacetime