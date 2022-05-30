/**
 * @file rotating-get_staples.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief
 * @version 0.1
 * @date 2022-05-30
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include "accum_type.hh"
#include "adjointfield.hh"
#include "gaugeconfig.hh"
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
   * @brief sum of all the staples, including the effect of metric
   * sum of regular and chair staples, weighted with their spacetime weights from the
   * metric g_{\mu\nu}
   * @tparam Group
   * @param U
   * @param Omega
   * @return Group
   */
  template <class T, class G>
  T get_staples_with_st_fact(const gaugeconfig<G> &U,
                        const nd_max_arr<size_t> &x,
                        const size_t &mu,
                        const double &Omega) {
    T St1;
    for (size_t nu = 0; nu < U.getndims(); nu++) {
      if (mu == nu) {
        continue;
      }

      // "/4.0" --> each plaquette comes from an average over the clover
      const T S_tt = get_staple<T, G>(U, x, mu, nu, true, true) / 4.0;
      const T S_ff = get_staple<T, G>(U, x, mu, nu, false, false) / 4.0;

      // 2 contributions from the clover plaquette at x
      St1 += plaq_factor(x, mu, nu, Omega) * (S_tt + S_ff);

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

    T St2;
    for (size_t nu = mu + 1; nu < U.getndims() - 1; nu++) {
      for (size_t rho = 1; rho < U.getndims(); rho++) {
        if (rho == mu || rho == nu || rho == 3) {
          continue;
        }

        // "/8.0" --> each chair loop is divided by 8.0
        const T S_tt_t = get_chair_staple<T, G>(U, x, mu, nu, rho, {1, 1}, 1) / 8.0;
        const T S_tt_f = get_chair_staple<T, G>(U, x, mu, nu, rho, {1, 1}, 0) / 8.0;
        const T S_tf_t = get_chair_staple<T, G>(U, x, mu, nu, rho, {1, 0}, 1) / 8.0;
        const T S_tf_f = get_chair_staple<T, G>(U, x, mu, nu, rho, {1, 0}, 0) / 8.0;
        const T S_ft_t = get_chair_staple<T, G>(U, x, mu, nu, rho, {0, 1}, 1) / 8.0;
        const T S_ft_f = get_chair_staple<T, G>(U, x, mu, nu, rho, {0, 1}, 0) / 8.0;
        const T S_ff_t = get_chair_staple<T, G>(U, x, mu, nu, rho, {0, 0}, 1) / 8.0;
        const T S_ff_f = get_chair_staple<T, G>(U, x, mu, nu, rho, {0, 0}, 0) / 8.0;

        // 4 contributions from the chairs at x
        St2 +=
          chair_stm_factor(x, mu, nu, rho, Omega) * (S_tt_t + S_ft_f - S_tf_f - S_ff_t);

        // x+\nu
        St2 += chair_stm_factor(xp(x, nu), mu, nu, rho, Omega) * (S_tt_t - S_tf_f);
        // x+\mu+\nu
        St2 +=
          chair_stm_factor(xp(xp(x, mu), nu), mu, nu, rho, Omega) * (S_tf_t - S_tt_f);
        // x+\mu
        St2 += chair_stm_factor(xp(x, mu), mu, nu, rho, Omega) *
               (S_tf_t + S_ff_f - S_tt_f - S_ft_t);
        // // x-\nu
        St2 += chair_stm_factor(xm(x, nu), mu, nu, rho, Omega) * (S_ft_f - S_ff_t);
        // x+\mu-\nu
        St2 +=
          chair_stm_factor(xp(xm(x, nu), mu), mu, nu, rho, Omega) * (S_ff_f - S_ft_t);
      }
    }

    T Stap = U(x, mu) * (St1 + St2);
    double bNc = U.getBeta() / double(U.getNc());
    return bNc * Stap;
  }

} // namespace rotating_spacetime
