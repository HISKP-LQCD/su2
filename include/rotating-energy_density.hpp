/**
 * @file rotating-energy_density.hpp
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief energy density in a rotating frame of reference
 * @version 0.1
 * @date 2022-05-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include "accum_type.hh"
#include "gaugeconfig.hh"
#include "tensors.hh"

#include "rotating-gauge_energy.hpp"

#include <math.h>

namespace rotating_spacetime {
  template <class T> using nd_max_arr = spacetime_lattice::nd_max_arr<T>;

  /**
   * @brief updating energy density and topological charge
   * this function is the generalization, to rotating frames of references, of the one in
   * flat spacetime
   * @tparam T
   * @param U
   * @param Omega
   * @param res
   * @param Q
   * @param cloverdef
   */
  template <class T>
  void energy_density(const gaugeconfig<T> &U,
                      const double &Omega,
                      double &res,
                      double &Q,
                      bool cloverdef = true) {
    res = 0.;
    Q = 0.;

    typedef typename accum_type<T>::type accum;
    // Euclidean 4D totally anti-symemtric tensor
    static epsilon4_t eps4 = new_epsilon4();

    nd_max_arr<size_t> x = {0, 0, 0, 0};
    for (x[0] = 0; x[0] < U.getLt(); x[0]++) {
      for (x[1] = 0; x[1] < U.getLx(); x[1]++) {
        for (x[2] = 0; x[2] < U.getLy(); x[2]++) {
          for (x[3] = 0; x[3] < U.getLz(); x[3]++) {
            nd_max_arr<size_t> x1 = x;
            nd_max_arr<size_t> x2 = x;
            nd_max_arr<size_t> x3 = x;
            accum G[4][4];
            for (size_t mu = 0; mu < U.getndims() - 1; mu++) {
              for (size_t nu = mu + 1; nu < U.getndims(); nu++) {
                accum leaf =
                  rotating_spacetime::plaquette<accum, T>(U, x, mu, nu, true, true);

                if (cloverdef) {
                  leaf = rotating_spacetime::clover_leaf<accum,T, nd_max_arr<size_t>>(U, x, mu, nu);
                }

                // traceless and anti-hermitian
                // here we include a factor 1/2 already
                G[mu][nu] = traceless_antiherm(leaf);
                // trace(G_{mu,nu}^a G_{mu,nu}^a)
                res += retrace(G[mu][nu] * G[mu][nu]);
              }
            }

            if (U.getndims() == 4) {
              // sum up the topological charge contribution now
              for (int i = 0; i < eps4.N; i++) {
                int i1 = eps4.eps_idx[i][0];
                int i2 = eps4.eps_idx[i][1];
                int i3 = eps4.eps_idx[i][2];
                int i4 = eps4.eps_idx[i][3];

                // when Gmunu components from the lower triangle are to be used,
                // we can simply skip them and multiply our normalisation by a factor of
                // four in total
                if (i2 < i1) {
                  continue;
                }
                if (i4 < i3) {
                  continue;
                }
                Q += eps4.eps_val[i] * retrace(G[i1][i2] * G[i3][i4]);
              }
            }
            if (U.getndims() == 2) {
              Q += -std::imag(trace((G[0][1] - G[1][0])));
            }
          }
        }
      }
    }
    // now we need to divide by 2, but we get a factor of two since we only
    // averaged mu < nu
    res = -res / U.getVolume();
    // averaged over four plaquette Wilson loops 1./4./4.
    // if (cloverdef)
    //   res /= 16.;
    // factor 4 from summing only mu < nu and rho < sigma
    Q = -4. * Q / (32.0 * M_PI * M_PI);
    // factor 1/16 from G_\mu\nu with clover definition
    if (cloverdef)
      Q /= 16.;
  }

} // namespace rotating_spacetime