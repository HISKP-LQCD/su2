/**
 * @file obc_weights.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief routines for imposing open boundary conditions on the gauge configuration
 * @version 0.1
 * @date 2023-01-30
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "gaugeconfig.hh"
#include "geometry.hh"

namespace obc {

  template <class Group> void apply_spatial_obc(gaugeconfig<Group> &U) {
    const size_t ndims = U.getndims();
    for (size_t x0 = 0; x0 < U.getLt(); x0++) {
      for (size_t x1 = 0; x1 < U.getLx(); x1++) {
        for (size_t x2 = 0; x2 < U.getLy(); x2++) {
          for (size_t x3 = 0; x3 < U.getLz(); x3++) {
            const std::vector<size_t> x = {x0, x1, x2, x3};
            const bool b1 = (x1 == 0 && ndims > 1);
            const bool b2 = (x2 == 0 && ndims > 2);
            const bool b3 = (x3 == 0 && ndims > 3);
            for (size_t mu = 0; mu < U.getndims(); mu++) {
              if (b1 || b2 || b3) {
                U(x, mu) = 0.0;
              }
            }
          }
        }
      }
    }
  }

} // namespace obc
