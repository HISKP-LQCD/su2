// Copyright (C) 2022 S. Romiti

#pragma once

#include "accum_type.hh"
#include "gaugeconfig.hh"
#include "geometry.hh"
#include "su2.hh"

#ifdef _USE_OMP_
#include <omp.h>
#endif
#include <fstream>
#include <iomanip>
#include <vector>

/**
 * @brief Polyakov loop at spatial point \vec{x}: P(\vec{x})
 * (see eq. 3.60 of https://link.springer.com/book/10.1007/978-3-642-01850-3)
 *
 * @tparam Group
 * @param U gauge configuration pointer
 * @param xi spatial components of the position
 * @return double
 */
template <class Group = su2>
double polyakov_loop(const gaugeconfig<Group> &U,
                     const std::array<int, spacetime_lattice::nd_max - 1> &xi) {
  typedef typename accum_type<Group>::type accum;

  const size_t mu = 0; // Polyakov loop contains U_{\mu=0}(x) only
  std::vector<size_t> x = {0, xi[0], xi[1], xi[2]};
  accum P(1., 0.);
  for (x[0] = 0; x[0] < U.getLt(); x[0]++) {
      P *= U(x, mu);
  }
  return retrace(P); // taking the real part averages over the 2 orientations
}

/**
 * @brief spatial average of the Polyakov loop
 * (see eq. 11 of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.103.094515)
 * @tparam Group
 * @param U
 * @return double
 */
template <class Group = su2> double polyakov_loop_spatial_average(gaugeconfig<Group> &U) {
  double ploop = 0.;
  typedef typename accum_type<Group>::type accum;
  const size_t ndims = U.getndims();

  const size_t nd_s = spacetime_lattice::nd_max - 1;
#pragma omp parallel for reduction(+: ploop)
  for (size_t x1 = 0; x1 < U.getLx(); x1++) {
    for (size_t x2 = 0; x2 < U.getLy(); x2++) {
      for (size_t x3 = 0; x3 < U.getLz(); x3++) {
        std::array<int, nd_s> x = {x1, x2, x3}; // spatial position
        ploop += polyakov_loop(U, x);
      }
    }
  }
  return ploop / U.getVolume();
}
