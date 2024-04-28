/**
 * @file uniform_sweeps.hpp
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief Sweeps for nested sampling
 *
 */

#pragma once

#include "accum_type.hh"
#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "random_element.hh"

// #ifdef _USE_OMP_
// #include <omp.h>
// #endif

#include <random>
#include <vector>

// apply n_sweeps sweeps drawing elements uniformly with the constraint of P < Pmax
template <class URNG, class Group>
void uniform_sweeps(gaugeconfig<Group> &U,
                    const double &P0,
                    const double &Pmax,
                    URNG engine,
                    const double &delta,
                    const size_t &n_sweeps) {
  std::uniform_real_distribution<double> uniform(0., 1.);
  typedef typename accum_type<Group>::type accum;

  const double norm_fact =
    double(U.getVolume() * U.getNc() * spacetime_lattice::num_pLloops_half(U.getndims()));

  double P_in = P0; // initial value of Plaquette average

  size_t i_sweep = 0;
  for (size_t x0 = 0; x0 < U.getLt(); x0++) {
    for (size_t x1 = 0; x1 < U.getLx(); x1++) {
      for (size_t x2 = 0; x2 < U.getLy(); x2++) {
        for (size_t x3 = 0; x3 < U.getLz(); x3++) {
          std::vector<size_t> x = {x0, x1, x2, x3};
          for (size_t mu = 0; mu < U.getndims(); mu++) {
            if (i_sweep > n_sweeps) {
              break;
            }
            Group R;
            accum K;
            get_staples_MCMC_step(K, U, x, mu, 1.0, false);
            random_element(R, engine, delta);
            double deltaP = retrace(U(x, mu) * R * K) - retrace(U(x, mu) * K);
            double P_new = P_in + (deltaP / norm_fact);
            if (P_new < Pmax) {
              U(x, mu) = U(x, mu) * R;
              P_in = P_new;
            }
            ++i_sweep;
          }
        }
      }
    }
  }
  return;
}
