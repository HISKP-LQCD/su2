/**
 * @file rotating-sweep.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief
 * @version 0.1
 * @date 2022-05-30
 *
 * @copyright Copyright (c) 2022
 */

#pragma once

#include "accum_type.hh"
#include "gaugeconfig.hh"
#include "random_element.hh"
#include "rotating-get_staples.hpp"

#ifdef _USE_OMP_
#include <omp.h>
#endif

#include <random>
#include <vector>

namespace rotating_spacetime {
  template <class T> using nd_max_arr = spacetime_lattice::nd_max_arr<T>;

  /**
   * @brief N_hit Metropolis-Updates
   * analogous of flat spacetime version, but in presence of a rotational metric
   * @param Omega (imaginary) angular frequency
   */
  template <class URNG, class Group>
  std::vector<double> sweep(gaugeconfig<Group> &U,
                            const double &Omega,
                            std::vector<URNG> engine,
                            const double &delta,
                            const size_t &N_hit,
                            const double &beta,
                            const double &xi = 1.0,
                            const bool &anisotropic = false) {
    std::uniform_real_distribution<double> uniform(0., 1.);
    typedef typename accum_type<Group>::type accum;
    size_t rate = 0, rate_time = 0;
#pragma omp parallel
    {
#ifdef _USE_OMP_
      size_t thread_num = omp_get_thread_num();
#else
      size_t thread_num = 0;
#endif
      for (size_t x0_start = 0; x0_start < 1; x0_start++) {
#pragma omp for reduction(+ : rate, rate_time)
        for (size_t x0 = 0; x0 < U.getLt(); x0 += 2) {
          Group R;
          for (size_t x1 = 0; x1 < U.getLx(); x1++) {
            for (size_t x2 = 0; x2 < U.getLy(); x2++) {
              for (size_t x3 = 0; x3 < U.getLz(); x3++) {
                nd_max_arr<size_t> x = {x0, x1, x2, x3};
                for (size_t mu = 0; mu < U.getndims(); mu++) {
                  accum K = rotating_spacetime::get_staples_with_st_fact<accum, Group>(U, x, mu, Omega);
                  for (size_t n = 0; n < N_hit; n++) {
                    random_element(R, engine[thread_num], delta);
                    double deltaS = beta / static_cast<double>(U.getNc()) *
                                    (retrace(U(x, mu) * K) - retrace(U(x, mu) * R * K));
                    bool accept = (deltaS < 0);
                    if (!accept)
                      accept = (uniform(engine[thread_num]) < exp(-deltaS));
                    if (accept) {
                      U(x, mu) = U(x, mu) * R;
                      U(x, mu).restoreSU();
                      rate += 1;
                      if (mu == 0) {
                        rate_time += 1;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    std::vector<double> res = {double(rate) / double(N_hit) / double(U.getSize()),
                               double(rate_time) / double(N_hit) / double(U.getVolume())};
    return res;
  }


} // namespace rotating_spacetime
