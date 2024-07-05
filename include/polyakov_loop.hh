// Copyright (C) 2022 S. Romiti

#pragma once
#include <complex>
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

using Complex = std::complex<double>;

namespace polyakov_loop{
/**
 * @brief Polyakov loop at spatial point \vec{x}: P(\vec{x})
 * (see eq. 3.60 of https://link.springer.com/book/10.1007/978-3-642-01850-3)
 *
 * @tparam Group
 * @param U gauge configuration pointer
 * @param xi spatial components of the position
 * @return Complex
 */
template <class Group = su2>
Complex polyakov_loop(const gaugeconfig<Group> &U,
                     const std::array<size_t, spacetime_lattice::nd_max - 1> &xi) {
  typedef typename accum_type<Group>::type accum;

  const size_t mu = 0; // Polyakov loop contains U_{\mu=0}(x) only
  std::vector<size_t> x = {0, xi[0], xi[1], xi[2]};
  //accum P;
  Group P;
  P.set_to_identity();

  for (x[0] = 0; x[0] < U.getLt(); x[0]++) {
      P *= U(x, mu);
  }
  return trace(P); // taking the rtrace part averages over the 2 orientations
}

/**
 * @brief Implements the correlator 3.61 in Gattringer Lang by fixing x at zero
*and averaging over x + r in all space directions
 * @tparam Group  
 * @param U configuration pointer
 * @param r Distance between x and y
 * @return Complex 
 */
template <class Group = su2>
Complex polyakov_loop_correlator(const gaugeconfig<Group> &U,
                                const size_t &r){
                                  const size_t nd_s = spacetime_lattice::nd_max - 1;
                                  std::array <size_t, nd_s> x = {0,0,0}; // x fixed at zero
                                  Complex polyakov_loop_x = conj(polyakov_loop(U, x)); //complex conjugate of polyakov loop for x 
                                  std::array <size_t, nd_s> y = {r, 0,0}; // x +r in x direction
                                  Complex correlator = 0.;

                                  correlator += polyakov_loop_x*polyakov_loop(U, y);
                                  y = {0,r,0}; //x + r in y direction
                                  correlator += polyakov_loop_x*polyakov_loop(U, y);
                                  y = {0,0,r}; // x + r in z direction
                                  correlator += polyakov_loop_x * polyakov_loop(U,y);
                                  correlator /= 3.; // average over all directions
                                  return correlator;
                                }

/**
 * @brief spatial average of the Polyakov loop
 * (see eq. 11 of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.103.094515)
 * @tparam Group
 * @param U
 * @return Complex
 */
template <class Group = su2>
 Complex polyakov_loop_spatial_average(const gaugeconfig<Group> &U) {
  Complex ploop = 0.;
  //typedef typename accum_type<Group>::type accum;
  const size_t ndims = U.getndims();

  const size_t nd_s = spacetime_lattice::nd_max - 1;
#pragma omp parallel for reduction(+: ploop)
  for (size_t x1 = 0; x1 < U.getLx(); x1++) {
    for (size_t x2 = 0; x2 < U.getLy(); x2++) {
      for (size_t x3 = 0; x3 < U.getLz(); x3++) {
        std::array<size_t, nd_s> x = {x1, x2, x3}; // spatial position
        ploop += polyakov_loop(U, x);
      }
    }
  }

  Complex helpvar = U.getVolume();
  return ploop / helpvar;
}
}