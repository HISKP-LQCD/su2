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
    int flag = 1;
    if (wx == 0.)
      flag = 0;
    for (size_t direction = 0; direction < r.size(); direction++) {
      directionloop = direction;//(direction + U.getndims()) % U.getndims();
      if (wx == 0.)
        flag = 0;
      for (size_t length = 0; length < r[direction]; length++) {
        if (wx == 0.)
          flag = 0;
        L *= U(xrun, directionloop);
        xrun[directionloop] += 1;
        wx = w(xrun);
      }
    }
    if (wx == 0.)
      flag = 0;
    for (size_t direction = 0; direction < r.size(); direction++) {
      directionloop = direction;//(direction + U.getndims()) % U.getndims();
      if (wx == 0.)
        flag = 0;
      for (size_t length = 0; length < r[direction]; length++) {
        xrun[directionloop] -= 1;
        wx = w(xrun);
        if (wx == 0.)
          flag = 0;
        L *= U(xrun, directionloop).dagger();
      }
    }
    double loop = double(flag) * retrace(L);
    return loop;
  }

} // namespace obc
