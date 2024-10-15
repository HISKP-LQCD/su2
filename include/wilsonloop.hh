// Copyright (C) 2017 C. Urbach

#pragma once

#include "accum_type.hh"
#include "gaugeconfig.hh"
#include "geometry.hh"
#include "su2.hh"
#include "partitionings.hh"

#ifdef _USE_OMP_
#include <omp.h>
#endif
#include <fstream>
#include <iomanip>
#include <vector>

/**
 * @brief Planar Wilson loop
 * Evaluation of the sum of all planar Wilson loop,
 * averaged over both orientations,
 * with the 1st vertex in any point of the lattice.
 * e.g. , starting from the origin:
 * (0,0) ->
 * (0, t*\hat{\nu}) ->
 * (t*\hat{\nu}, r*\hat{\mu} + t*\hat{nu}) ->
 * (r*\hat{\mu} + t*\hat{nu}, r*\hat{\mu}) ->
 * (r*\hat{\mu}, 0) -> (0,0)
 *
 * 's' ant 't' are meant to be respectively the spatial and temporal number of lattice
 * points
 *
 * Notes:
 *
 * 1. By traslational invariance, they're all the same analytically)
 *    See eq. (3.50) of https://link.springer.com/book/10.1007/978-3-642-01850-3 for
 * reference.
 * 2. Complex conjugation of links simply reverses the loop orientation
 *    [see e.g. eq. (2.34) and sec. 3.3.3 of
 * https://link.springer.com/book/10.1007/978-3-642-01850-3] This corresponds to the
 * interchange of quark and antiquark. Both orientations are valid in order to determine
 * the static potential.
 *
 * @tparam Group
 * @param U gauge configuration pointer
 * @param r number of steps in the \mu direction
 * @param t number of steps in the \nu direction
 * @param mu 1st direction of the loop
 * @param nu 2nd direction of the loop
 * @return double
 */
template <class Group = su2>
double planar_wilsonloop_dir(const gaugeconfig<Group> &U,
                             const size_t r,
                             const size_t t,
                             const size_t mu,
                             const size_t nu) {
  double loop = 0.;
  typedef typename accum_type<Group>::type accum;

  std::vector<size_t> x = {0, 0, 0, 0};
  for (x[0] = 0; x[0] < U.getLt(); x[0]++) {
    for (x[1] = 0; x[1] < U.getLx(); x[1]++) {
      for (x[2] = 0; x[2] < U.getLy(); x[2]++) {
        for (x[3] = 0; x[3] < U.getLz(); x[3]++) {
          std::vector<size_t> xrun = x;
          Group L;
          L.set_to_identity(); // L = 1.0
          for (size_t _t = 0; _t < t; _t++) {
            L *= U(xrun, nu);
            xrun[nu] += 1;
          }
          for (size_t s = 0; s < r; s++) {
            L *= U(xrun, mu);
            xrun[mu] += 1;
          }
          for (size_t _t = 0; _t < t; _t++) {
            xrun[nu] -= 1;
            L *= U(xrun, nu).dagger();
          }
          for (size_t s = 0; s < r; s++) {
            xrun[mu] -= 1;
            L *= U(xrun, mu).dagger();
          }
          loop += retrace(L); // taking the real part averages over the 2 orientations
        }
      }
    }
  }
  return loop;
}

/**
 * calculates the Wilson-loop given by the path in r
 * r[0] steps are taken in direction 0, r[1] steps in direction 1 and so on, with r[n]
 * steps taken in direction n%ndims For each direction for each step, the corresponding
 * link is multiplied onto the loop (standard Wilson-Loop definition): loop *=
 * prod_{i=0}^{r[n]} U_{n%ndims} (x+i*e_{n%ndims}+shifts from eaarlier steps) If the path
 * is done, it is traced back in the same direction, this time using the daggered links
 * The loop is calculated for each lattice point and averaged over the entire lattice
 * parallelization trivial
 * r=(1,1,0,0)=r(1,1) is the temporal plaquette, calculated in the order t->x.
 * the order x->t can be achieved by using r=(0,1,0,0,1) and 4d or r=(0,1,0,1) in 3d.
 * */
template <class Group = su2>
double wilsonloop_non_planar(const gaugeconfig<Group> &U, std::vector<size_t> r) {
  // goes path outlined in r in direction t->x->y->z, could go with other orders by using
  // longer vector r and inserting zeros parallelized with code from gauge_energy
  double loop = 0.;
  typedef typename accum_type<Group>::type accum;

#pragma omp parallel for reduction(+ : loop)
  for (size_t x0 = 0; x0 < U.getLt(); x0++) {
    for (size_t x1 = 0; x1 < U.getLx(); x1++) {
      for (size_t x2 = 0; x2 < U.getLy(); x2++) {
        for (size_t x3 = 0; x3 < U.getLz(); x3++) {
          std::vector<size_t> xrun = {x0, x1, x2, x3};
          Group L;
          L.set_to_identity(); // L = 1.0
          // needed if vector with directions contains more than 4 entries/if another
          // order than t-x-y-z is wanted
          size_t directionloop;
          for (size_t direction = 0; direction < r.size(); direction++) {
            directionloop = (direction + U.getndims()) % U.getndims();
            for (size_t length = 0; length < r[direction]; length++) {
              L *= U(xrun, directionloop);
              xrun[directionloop] += 1;
            }
          }
          for (size_t direction = 0; direction < r.size(); direction++) {
            directionloop = (direction + U.getndims()) % U.getndims();
            for (size_t length = 0; length < r[direction]; length++) {
              xrun[directionloop] -= 1;
              L *= U(xrun, directionloop).dagger();
            }
          }
          loop += retrace(L);
        }
      }
    }
  }
  return loop;
}

/**
 * @brief average of Wilson loops over all spatial directions
 * (0 is the temporal direction)
 * @tparam Group
 * @param U gauge configuration
 * @param r spatial extent of the loops
 * @param t temporal extent of the loops
 * @return double
 */
template <class Group>
double wilsonloop(const gaugeconfig<Group> &U, const size_t r, const size_t t) {
  double loop = 0.;
  const size_t ndims = U.getndims();
#pragma omp parallel for reduction(+ : loop)
  for (size_t mu = 1; mu < ndims; mu++) {
    loop += planar_wilsonloop_dir(U, r, t, mu, 0);
  }
  return loop / U.getVolume() / double(U.getNc()) / ndims;
}

/**
 * @brief saving all planar loop of the lattice grid
 * Computing and printing averages of all planar loops
 * for all possible spatial and temporal extents. It is assumed that Lx==Ly==Lz
 * @tparam Group
 * @param U gauge config
 * @param path path of the output file
 */
template <class Group>
void compute_all_loops(const gaugeconfig<Group> &U, std::string const &path) {
  // checking the spatial symmetry of the lattice
  const size_t Lt = U.getLt(), Lx = U.getLx();

  std::ostringstream oss;

  // printing a header specifying what columns contain
  oss << "t";
  for (size_t r = 1; r < Lx; r++) {
    oss << " r=" << r;
  }
  oss << "\n";

  // printing the data in the format t L(r=1) L(r=2) ... L(r=Lx-1)
  for (size_t t = 1; t < Lt; t++) {
    oss << t;
    for (size_t r = 1; r < Lx; r++) {
      double loop = wilsonloop(U, r, t);
      oss << " " << std::scientific << std::setprecision(15) << loop;
    }
    oss << std::endl;
  }

  std::ofstream ofs(path, std::ios::out);
  ofs << oss.str();

  return;
}

template <class Group>
void compute_spacial_loops(gaugeconfig<Group> &U, std::string const &path) {
  std::ofstream os(path, std::ios::out);
  size_t r[2] = {2, 8};
  for (size_t t = 1; t < U.getLt(); t++) {
    os << t;
    for (size_t i = 0; i < 2; i++) {
      double loop = wilsonloop(U, r[i], t);
      os << " " << std::scientific << std::setw(15) << loop;
    }
    os << std::endl;
  }
  return;
}
