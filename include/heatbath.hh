/**
 * @file heatbath.hpp
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
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
#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "random_element.hh"

#ifdef _USE_OMP_
#include <omp.h>
#endif

#include <random>
#include <vector>

namespace flat_spacetime {

  /**
   * @brief heatbath step
   *
   * Sampling each U_\mu(x) according to eq. (10) of
   * https://www.sciencedirect.com/science/article/pii/S0010465509002148   *
   *
   * @tparam URNG
   * @tparam Group
   * @param U gauge configuration
   * @param engine engine for RNG: size must be 2*number of threads
   * @param N_hit number of hits for the sampling
   * @param beta coupling beta in the action
   * @param xi bare anisotropy
   * @param anisotropic bool. flag, true when considering an anisotropic lattice
   * @return std::vector<double>: acceptance rates: {overall, only temporal ones}
   */
  template <class URNG, class Group>
  void heatbath(gaugeconfig<Group> &U,
                std::vector<URNG> engine,
                const size_t &N_hit,
                const double &beta,
                const double &xi = 1.0,
                const bool &anisotropic = false);

  template <class URNG>
  void heatbath(gaugeconfig<u1> &U,
                std::vector<URNG> engine,
                const size_t &N_hit,
                const double &beta,
                const double &xi = 1.0,
                const bool &anisotropic = false) {
    const double coupl_fact = (beta / double(U.getNc()));
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    typedef typename accum_type<u1>::type accum;

#ifdef _USE_OMP_
#pragma omp parallel
    {
      size_t thread_num = omp_get_thread_num();
#else
    size_t thread_num = 0;
#endif

      for (size_t x0_start = 0; x0_start < 2; x0_start++) {
#pragma omp for
        for (size_t x0 = x0_start; x0 < U.getLt(); x0 += 2) {
          for (size_t x1 = 0; x1 < U.getLx(); x1++) {
            for (size_t x2 = 0; x2 < U.getLy(); x2++) {
              for (size_t x3 = 0; x3 < U.getLz(); x3++) {
                const std::vector<size_t> x = {x0, x1, x2, x3};
                for (size_t mu = 0; mu < U.getndims(); mu++) {
                  accum K;
                  get_staples_MCMC_step(K, U, x, mu, xi, anisotropic);
                  const double theta_stap = get_phase(K);
                  const double rho = coupl_fact * get_abs(K);
                  const double A = rho / (2.0 * sinh(rho));
                  const double u1 = uniform(engine[thread_num]);
                  const double y = (1 / rho) * log((rho / A) * u1 + exp(-rho));
                  double phi = acos(y);
                  const double u2 = uniform(engine[thread_num]);
                  if (u2 >= 0.5) {
                    phi *= -1.0;
                  }
                  U(x, mu).set(phi - theta_stap);
                }
              }
            }
          }
        }
      }
#ifdef _USE_OMP_
    }
#endif
    return;
  }

  template <class URNG>
  void heatbath(gaugeconfig<su2> &U,
                std::vector<URNG> engine,
                const size_t &N_hit,
                const double &beta,
                const double &xi = 1.0,
                const bool &anisotropic = false) {
    spacetime_lattice::fatal_error("heatbath not implemented for SU(2)!", __func__);

    return;
  }

  template <class URNG>
  void heatbath(gaugeconfig<su3> &U,
                std::vector<URNG> engine,
                const size_t &N_hit,
                const double &beta,
                const double &xi = 1.0,
                const bool &anisotropic = false) {
    spacetime_lattice::fatal_error("heatbath not implemented for SU(3)!", __func__);

    return;
  }

} // namespace flat_spacetime