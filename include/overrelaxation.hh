/**
 * @file overrelaxation.hpp
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-05-30
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include "accum_type.hh"
#include "errors.hpp"
#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "random_element.hh"

#ifdef _USE_OMP_
#include <omp.h>
#endif

#include <random>
#include <vector>

/**
 * @brief  eq. below (4.50) of https://link.springer.com/book/10.1007/978-3-642-01850-3
 *
 * NOTE: the `engines` arguments is not used for U(1),
 * but is needed to have a consistent overflow with the other SU(N) groups
 */
template<class URNG>
void overrelaxation(gaugeconfig<u1> &U,
                    std::vector<URNG> engines,
                    const double &xi,
                    const bool &anisotropic) {
  typedef typename accum_type<u1>::type accum;

  const size_t endmu = U.getndims();
  for (size_t x0_start = 0; x0_start < 2; x0_start++) {
#pragma omp parallel for
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

/**
 * @brief overrelaxation step
 *
 * See eq. 7.34 of https://www.worldscientific.com/worldscibooks/10.1142/6065 for SU(N)
 */
template<class URNG>
void overrelaxation(gaugeconfig<su2> &U,
                    std::vector<URNG> engines,
                    const double &xi,
                    const bool &anisotropic) {
  typedef typename accum_type<su2>::type accum;

  const size_t endmu = U.getndims();
  for (size_t x0_start = 0; x0_start < 2; x0_start++) {
#pragma omp parallel for
    for (size_t x0 = x0_start; x0 < U.getLt(); x0 += 2) {
      size_t thread_num = omp_get_thread_num();
      for (size_t x1 = 0; x1 < U.getLx(); x1++) {
        for (size_t x2 = 0; x2 < U.getLy(); x2++) {
          for (size_t x3 = 0; x3 < U.getLz(); x3++) {
            const std::vector<size_t> x = {x0, x1, x2, x3};
            for (size_t mu = 0; mu < endmu; mu++) {
              su2 V;
              accum S;
              get_staples_MCMC_step(S, U, x, mu, xi, anisotropic);
              // NOTE: for SU(2) matrices the determinant is real
              const double detS = S.det().real();
              if (detS == 0.0) {
                // random element in the group
                random_element(V, engines[thread_num], 1.0);
              } else {
                // for SU(2) it is sufficient to generate V by normalizing the sum of
                // staples
                const double sqrt_det_S = sqrt(detS);
                S = (1.0 / sqrt(S.det())) * S;
                V.set(S.geta(), S.getb());
                V = V.dagger();
              }
              U(x, mu) = V * U(x, mu) * V;
            }
          }
        }
      }
    }
  }

  return;
}

template<class URNG>
void overrelaxation(gaugeconfig<su3> &U,
                    std::vector<URNG> engines,
                    const double &xi,
                    const bool &anisotropic) {
  fatal_error("overrelaxation not implemented for SU(3)!", __func__);

  return;
}
