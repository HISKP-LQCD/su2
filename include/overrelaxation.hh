/**
 * @file overrelaxation.hpp
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

/**
 * @brief overrelaxation step
 *
 * See eq. 7.34 of https://www.worldscientific.com/worldscibooks/10.1142/6065 for SU(N)
 *
 * The change in the action \Delta S accepted with probability min(1, exp(-\Delta S)).
 * @tparam URNG
 * @tparam Group
 * @param U gauge configuration
 * @param engine engine for random number generation
 * @param beta coupling beta in the action
 * @param xi bare anisotropy
 * @param anisotropic bool. flag, true when considering an anisotropic lattice
 * @return std::vector<double>: acceptance rates: {overall, only temporal ones}
 */
template <class Group>
void overrelaxation(gaugeconfig<Group> &U,
                    const double &beta,
                    const double &xi = 1.0,
                    const bool &anisotropic = false,
                    const bool & temporalonly = false);

/**
 * @brief  eq. below (4.50) of https://link.springer.com/book/10.1007/978-3-642-01850-3
 */
template <>
void overrelaxation(gaugeconfig<u1> &U,
                    const double &beta,
                    const double &xi,
                    const bool &anisotropic,
                    const bool & temporalonly) {
  typedef typename accum_type<u1>::type accum;
  const size_t endmu = temporalonly ? 1 : U.getndims(); 

  for (size_t x0_start = 0; x0_start < 2; x0_start++) {
#pragma omp for
    for (size_t x0 = x0_start; x0 < U.getLt(); x0 += 2) {
      for (size_t x1 = 0; x1 < U.getLx(); x1++) {
        for (size_t x2 = 0; x2 < U.getLy(); x2++) {
          for (size_t x3 = 0; x3 < U.getLz(); x3++) {
            const std::vector<size_t> x = {x0, x1, x2, x3};
            for (size_t mu = 0; mu < endmu; mu++) {
              accum K;
              get_staples_MCMC_step(K, U, x, mu, xi, anisotropic);
              const double phi = get_phase(K);
              U(x, mu).set(-2 * phi - U(x, mu).geta());
            }
          }
        }
      }
    }
  }

  return;
}

template <>
void overrelaxation(gaugeconfig<su2> &U,
                    const double &beta,
                    const double &xi,
                    const bool &anisotropic,
                    const bool & temporalonly) {
  spacetime_lattice::fatal_error("overrelaxation not implemented for SU(2)!", __func__);

  return;
}

template <>
void overrelaxation(gaugeconfig<su3> &U,
                    const double &beta,
                    const double &xi,
                    const bool &anisotropic,
                    const bool & temporalonly) {
  spacetime_lattice::fatal_error("overrelaxation not implemented for SU(3)!", __func__);

  return;
}
