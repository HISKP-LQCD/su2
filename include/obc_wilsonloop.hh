// wilson loops with open boundary conditions

#pragma once

#include "accum_type.hh"
#include "gaugeconfig.hh"
#include "geometry.hh"
#include "obc_weights.hh"
#include "su2.hh"

#ifdef _USE_OMP_
#include <omp.h>
#endif
#include <fstream>
#include <iomanip>
#include <vector>

namespace obc {

  /**
   * @brief Planar Wilson loop
   *
   * The only difference with the homonymous routine are the boundary conditions.
   * See the latter documentation.
   */
  template <class Group = su2>
  double planar_wilsonloop_dir(const gaugeconfig<Group> &U,
                               const obc::weights &w,
                               const std::vector<size_t> &x,
                               const size_t &r,
                               const size_t &t,
                               const size_t &mu,
                               const size_t &nu) {
    typedef typename accum_type<Group>::type accum;

    std::vector<size_t> xrun = x;
    const double wx = w(xrun); // changes along the loop
    Group L;
    L.set_to_identity(); // L = 1.0
    for (size_t _t = 0; _t < t; _t++) {
      L *= wx * U(xrun, nu);
      xrun[nu] += 1;
      wx = w(xrun);
    }
    for (size_t s = 0; s < r; s++) {
      L *= U(xrun, mu);
      xrun[mu] += 1;
      wx = w(xrun);
    }
    for (size_t _t = 0; _t < t; _t++) {
      xrun[nu] -= 1;
      wx = w(xrun);
      L *= U(xrun, nu).dagger();
    }
    for (size_t s = 0; s < r; s++) {
      xrun[mu] -= 1;
      wx = w(xrun);
      L *= U(xrun, mu).dagger();
    }
    double loop = retrace(L); // taking the real part averages over the 2 orientations
    return loop;
  }

  /**
   * @brief Non planar Wilson loop
   *
   * The only difference with the homonymous routine are the boundary conditions.
   * See the latter documentation.
   *
   * */
  template <class Group = su2>
  double wilsonloop_non_planar(const gaugeconfig<Group> &U,
                               const obc::weights &w,
                               const std::vector<size_t> &x,
                               const std::vector<size_t> &r) {
    typedef typename accum_type<Group>::type accum;

    std::vector<size_t> xrun = x;
    Group L;
    L.set_to_identity(); // L = 1.0
    // needed if vector with directions contains more than 4 entries/if another
    // order than t-x-y-z is wanted
    size_t directionloop;
    double wx = w(xrun);
    for (size_t direction = 0; direction < r.size(); direction++) {
      directionloop = (direction + U.getndims()) % U.getndims();
      for (size_t length = 0; length < r[direction]; length++) {
        L *= U(xrun, directionloop);
        xrun[directionloop] += 1;
        wx = w(xrun);
      }
    }
    for (size_t direction = 0; direction < r.size(); direction++) {
      directionloop = (direction + U.getndims()) % U.getndims();
      for (size_t length = 0; length < r[direction]; length++) {
        xrun[directionloop] -= 1;
        wx = w(xrun);
        L *= U(xrun, directionloop).dagger();
      }
    }
    double loop = retrace(L);
    return loop;
  }

} // namespace obc
