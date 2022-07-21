/**
 * @file glueballs.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief routines for the calculation of the glueball correlators
 * @version 0.1
 * @date 2022-07-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include "operators.hpp"

namespace glueballs {

  template <class T> using nd_max_arr = spacetime_lattice::nd_max_arr<T>;

  /**
   * @brief Rectangular Wilson loop of size (a,b) in the (mu,nu) plane
   * Px : boolean flag for the application of the parity operator
   */
  template <class Float, class Group>
  std::complex<Float> trace_wloop_munu(const gaugeconfig<Group> &U,
                                       nd_max_arr<int> x,
                                       const size_t &a,
                                       const size_t &b,
                                       const size_t &mu,
                                       const size_t &nu,
                                       const bool &Px) {
    typedef typename accum_type<Group>::type accum;

    accum L = U(x, mu) * (U(x, mu).dagger()); // "1", independently of the group
    for (size_t s = 0; s < a; s++) {
      L *= operators::parity(Px, 0, U, x, mu);
      x[mu] += 1;
    }
    for (size_t _t = 0; _t < b; _t++) {
      L *= operators::parity(Px, 0, U, x, nu);
      x[nu] += 1;
    }
    for (size_t s = 0; s < a; s++) {
      x[mu] -= 1;
      L *= operators::parity(Px, 0, U, x, mu).dagger();
    }
    for (size_t _t = 0; _t < b; _t++) {
      x[nu] -= 1;
      L *= operators::parity(Px, 0, U, x, nu).dagger();
    }

    return trace(L);
  }

  /**
   * @brief trace_wloop at \vec{p}=\vec{0}
   */
  template <class Float, class Group>
  std::complex<Float> rest_trace_wloop_munu(size_t &t,
                                            const gaugeconfig<Group> &U,
                                            const size_t &a,
                                            const size_t &b,
                                            const size_t &mu,
                                            const size_t &nu,
                                            const bool &Px) {
    nd_max_arr<int> x = {int(t), 0, 0, 0};
    std::complex<Float> loop = 0.0;
    for (x[1] = 0; x[1] < U.getLx(); x[1]++) {
      for (x[2] = 0; x[2] < U.getLy(); x[2]++) {
        for (x[3] = 0; x[3] < U.getLz(); x[3]++) {
          loop += trace_wloop_munu<Float, Group>(U, x, a, b, mu, nu, Px);
        }
      }
    }
    const double spat_vol = double(U.getVolume()) / double(U.getLt());
    return loop / spat_vol;
  }

} // namespace glueballs
