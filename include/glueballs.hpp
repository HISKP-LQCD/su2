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
   * @brief L-shaped Wilson loop of size (a,b) in the (mu,nu) directions
   * @param a shortest size on the mu direction
   * @param b shortest size on the nu direction
   * @param Px : boolean flag for the application of the parity operator
   */
  template <class Float, class Group>
  std::complex<Float> trace_Lloop_munu(const gaugeconfig<Group> &U,
                                       nd_max_arr<int> x,
                                       const size_t &a,
                                       const size_t &b,
                                       const size_t &mu,
                                       const size_t &nu,
                                       const bool &Px) {
    typedef typename accum_type<Group>::type accum;

    accum L(1.0);
    // start at x
    // go to x + 2*a*\hat{mu}
    for (size_t s = 0; s < 2 * a; s++) {
      L *= operators::parity(Px, 0, U, x, mu);
      x[mu] += 1;
    }
    // go to x + 2*a*\hat{mu} + 2*b*\hat{nu}
    for (size_t _t = 0; _t < 2 * b; _t++) {
      L *= operators::parity(Px, 0, U, x, nu);
      x[nu] += 1;
    }
    // go to x + a*\hat{mu} + 2*b*\hat{nu}
    for (size_t s = 0; s < a; s++) {
      x[mu] -= 1;
      L *= operators::parity(Px, 0, U, x, mu).dagger();
    }
    // go to x + a*\hat{mu} + b*\hat{nu}
    for (size_t _t = 0; _t < b; _t++) {
      x[nu] -= 1;
      L *= operators::parity(Px, 0, U, x, nu).dagger();
    }
    // go to x + b*\hat{nu}
    for (size_t s = 0; s < a; s++) {
      x[mu] -= 1;
      L *= operators::parity(Px, 0, U, x, mu).dagger();
    }
    // go to x
    for (size_t _t = 0; _t < b; _t++) {
      x[nu] -= 1;
      L *= operators::parity(Px, 0, U, x, nu).dagger();
    }
    // loop closed

    return trace(L);
  }

  /**
   * @brief trace_wloop_munu at \vec{p}=\vec{0}
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

  /**
   * @brief trace_Lloop_munu at \vec{p}=\vec{0}
   */
  template <class Float, class Group>
  std::complex<Float> rest_trace_Lloop_munu(size_t &t,
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
          loop += trace_Lloop_munu<Float, Group>(U, x, a, b, mu, nu, Px);
        }
      }
    }
    const double spat_vol = double(U.getVolume()) / double(U.getLt());
    return loop / spat_vol;
  }

  /**
   * @brief average over dimensions of trace_wloop at \vec{p}=\vec{0}
   * average of rectangular wilson loops of size a \times b over all dimensions.
   * @param statial true when the sum is done only on the i,j loops
   */
  template <class Float, class Group>
  std::complex<Float> rest_trace_wloop(size_t &t,
                                       const gaugeconfig<Group> &U,
                                       const size_t &a,
                                       const size_t &b,
                                       const bool &Px,
                                       const bool &spatial = false) {
    const size_t ss = size_t(spatial);
    const size_t d = U.getndims();
    std::complex<Float> loop = 0.0;
    for (size_t mu = ss; mu < d; mu++) {
      for (size_t nu = mu + 1; nu < d; nu++) {
        loop += rest_trace_wloop_munu<Float, Group>(t, U, a, b, mu, nu, Px);
      }
    }
    const double norm_fact = spacetime_lattice::num_pLloops_half(d - ss);
    return loop / norm_fact;
  }

  /**
   * @brief average over dimensions of trace_Lloop at \vec{p}=\vec{0}
   * average of L-shaped wilson loops of size a \times b over all dimensions.
   * @param statial true when the sum is done only on the i,j loops
   */
  template <class Float, class Group>
  std::complex<Float> rest_trace_Lloop(size_t &t,
                                       const gaugeconfig<Group> &U,
                                       const size_t &a,
                                       const size_t &b,
                                       const bool &Px,
                                       const bool &spatial = false) {
    const size_t ss = size_t(spatial);
    const size_t d = U.getndims();
    std::complex<Float> loop = 0.0;
    for (size_t mu = ss; mu < d; mu++) {
      for (size_t nu = mu + 1; nu < d; nu++) {
        loop += rest_trace_Lloop_munu<Float, Group>(t, U, a, b, mu, nu, Px);
      }
    }
    const double norm_fact = spacetime_lattice::num_pLloops_half(d - ss);
    return loop / norm_fact;
  }

  template <class Float, class Group>
  std::complex<Float> get_rest_trace_loop(const std::string type,
                                          size_t &t,
                                          const gaugeconfig<Group> &U,
                                          const size_t &a,
                                          const size_t &b,
                                          const bool &Px,
                                          const bool &spatial = false) {
    if (type == "rectangle") {
      return rest_trace_wloop<Float, Group>(t, U, a, b, Px, spatial);
    } else if (type == "L") {
//      return rest_trace_Lloop<Float, Group>(t, U, a, b, Px, spatial);
      return rest_trace_Lloop_munu<Float, Group>(t, U, a, b, 1, 2, Px);
    } else {
      std::cerr << "Error: Invalid loop type '" << type << "'.\nAborting.\n";
      std::abort();
    }
  }

} // namespace glueballs
