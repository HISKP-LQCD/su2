#pragma once

#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "random_element.hh"
#include "su2.hh"
#include "u1.hh"
#include "su3.hh"
#include "su3_accum.hh"

#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#ifdef _USE_OMP_
#include <omp.h>
#endif

// prescription for smearing taken from https://arxiv.org/abs/hep-lat/0209159
// not tested for SU2, for U1 lattice stays the same when smearing with alpha=1
// norm is the number of staples which are used in the summation
// if spacial=true, only the staples with mu>=1 are smeared and taken into account for the
// staples
template <class Group>
void smearlatticeape(gaugeconfig<Group> &U,
                     double alpha,
                     bool spatial_only = false,
                     bool temporal_only = false) {
  if (U.getndims() == 2 && spatial_only) {
    std::cerr << "Spacial smearing is not possible in two dimensions!" << std::endl;
    return;
  }
  if (spatial_only && temporal_only) {
    std::cerr
      << "selecting spatial_only and temporal_only at the same time does not make sense!"
      << std::endl;
    return;
  }
  gaugeconfig<Group> Uold = U;
  typedef typename accum_type<Group>::type accum;
  size_t startmu = 0;
  size_t endmu = U.getndims();
  size_t norm = 2.0 * (U.getndims() - 1);
  if (spatial_only) {
    startmu = 1;
    norm = 2.0 * (U.getndims() - 2);
  }
  if (temporal_only) {
    startmu = 0;
    endmu = 1;
  }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (size_t x0 = 0; x0 < U.getLt(); x0++) {
    for (size_t x1 = 0; x1 < U.getLx(); x1++) {
      for (size_t x2 = 0; x2 < U.getLy(); x2++) {
        for (size_t x3 = 0; x3 < U.getLz(); x3++) {
          std::vector<size_t> x = {x0, x1, x2, x3};
          for (size_t mu = startmu; mu < endmu; mu++) {
            // K is intialized to (0,0) even if not explicitly specified
            accum K(0.0, 0.0);
            get_staples_MCMC_step(K, Uold, x, mu, 1.0, false, spatial_only);
            Group Uprime(Uold(x, mu) * alpha + K * (1 - alpha) / double(norm));
            U(x, mu) = Uprime;
            U(x, mu).restoreSU();
          }
        }
      }
    }
  }
}

/**
 * @brief wrapper for smearlatticeape with spatial_only==true, temporal_only==false
 */
template <class Group> void spatial_smearlatticeape(gaugeconfig<Group> &U, double alpha) {
  smearlatticeape(U, alpha, true, false);
}

/**
 * @brief smearing with prescription generalised from https://arxiv.org/pdf/hep-lat/0209159.pdf
 * default: smearing onlz of spatial links with spatial staples
 * number of staples in 'positive' direction: 
 * d 
 * -1 (direction of original link is not used for staple) 
 * -1 if spatial (temporal links not taken into account)
 * have to multiply by 2 to also account for 'negative' direction
 * beta(spatial) = (1-alpha)/(# staples) = (1-alpha)/(2*(d-1-1)) =(1-alpha)/(2*d-4)
 * prescription for smearing: U <- alpha*U + beta*sum(staples)
 */
template <class Float, class Group>
void APEsmearing(gaugeconfig<Group> &U, const double &alpha, const bool spatial=true) {
  size_t d = U.getndims();

  if (d == 2) {
    std::cerr << "Spatial smearing is not possible in 2 dimensions!" << std::endl;
    return;
  }
  const gaugeconfig<Group> Uold = U;
  typedef typename accum_type<Group>::type accum;
  size_t startmu = size_t(spatial); // 0 or 1, casted from bool
  const Float beta = (1.0-alpha) / (double(2*(d-1-startmu)));

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (size_t x0 = 0; x0 < U.getLt(); x0++) {
    for (size_t x1 = 0; x1 < U.getLx(); x1++) {
      for (size_t x2 = 0; x2 < U.getLy(); x2++) {
        for (size_t x3 = 0; x3 < U.getLz(); x3++) {
          std::vector<size_t> x = {x0, x1, x2, x3};
          for (size_t i = startmu; i < d; i++) {
            // K is intialized to (0,0) even if not explicitly specified
            accum K;
            get_staples_APE(K, Uold, x, i, spatial);
            K = alpha*accum(Uold(x, i)) + beta*K;
            const Group Uprime = accum_to_Group(K);
            U(x, i) = Uprime;
            U(x, i).restoreSU();
          }
        }
      }
    }
  }
  return;
}

// wrapper so the name of the other function can be changed without changing any other code
template <class Float, class Group>
void spatial_APEsmearing(gaugeconfig<Group> &U, const double &alpha) {
    APEsmearing<Float, Group>(U, alpha, true);
}

